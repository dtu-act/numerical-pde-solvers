classdef HelmholtzSolver2D < handle
    
    properties
        acoustic_params, simulation_params
        L, Ltilde, b, btilde, M, Sx, Sy
        P, x2D, y2D, conn, mesh_info
        boundary_type
        source
        mesh_created
        boundary_cond_f
    end

    methods
        function self = HelmholtzSolver2D(acoustic_params)
            self.acoustic_params = acoustic_params;
            self.mesh_created = false;
        end
        
        function setupMeshUniform(self, bbox, NeX, NeY)
            [XY, etov, bdetect_f] = meshing.mesh2D(bbox, NeX, NeY);
            self.mesh_info = sem.Mesh2D(XY,bbox,etov,bdetect_f);
            self.mesh_created = true;
        end

        function setupMeshNonUniform(self,bbox,xy0,hmin,hmax)
            [XY, etov, ~] = meshing.distMeshGeneration(bbox,xy0,hmin,hmax);          
            bound_detect_f = meshing.createBoundRectDetectFunc(bbox, 1e-10);            
            self.mesh_info = sem.Mesh2D(XY,bbox,etov,bound_detect_f);

            self.mesh_created = true;
        end

        function setupSolver(self, P_order, boundary_type, source, boundary_cond_f)            
            if ~self.mesh_created
                error('Mesh not created - call setupMeshUniform | setupMeshNonUniform')
            end
            self.clearResults()
            
            if nargin == 4
                self.boundary_cond_f = @(x,y) 0;
            else
                self.boundary_cond_f = boundary_cond_f;
            end

            self.boundary_type = boundary_type;
            self.source = source;            

            self.simulation_params = sem.SimulationParameters2D(P_order);
            Np = self.simulation_params.Np;
            Nfaces = self.simulation_params.Nfaces;
            
            etov = self.mesh_info.etov;
            mesh = self.mesh_info.XY;
            k = self.acoustic_params.k;

            Nk = size(etov,1); % number of elements (triangles)        

            vx = mesh(:,1);
            vy = mesh(:,2);

            [X,Y,r,s] = sem.verticesToNodes2D(P_order,etov,vx,vy);
            [X,Y,r,s] = sem.reordernodes2D(P_order,X,Y,r,s);

            % connectivity
            [etoe, etof] = spectral.tiConnect2D(etov, Nfaces);

            [self.conn, self.x2D, self.y2D, ~] = meshing.connectivityTable2D(P_order, X, Y, vx, vy, etov, etoe, etof);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            %%        ASSEMBLY        %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Build reference element matrices
            V = spectral.Vandermonde2D(P_order, r, s);
            [Dr, Ds] = spectral.Dmatrices2D(P_order, r, s, V);
            geofact = spectral.GeometricFactors2D(X, Y, Dr, Ds);

            % general mass matrix
            Me = inv(V*V');    

            [self.M, self.Sx, self.Sy] = sem.assembly.assemblyA2D(Nk, Np, geofact, self.conn, Me, Dr, Ds);

            switch source.type
                case models.SourceType.Function
                    % alternative: b_alt = Mxy*source.F(x2D(:),y2D(:));
                    self.b = sem.assembly.assemblyb2D(Nk, Np, geofact, self.conn, Me, self.x2D, self.y2D, source.F);
                case models.SourceType.PointSource
                    % F = dirac2D(source.r0(1),source.r0(2),source.Q,1e-10);
                    % alternative: b = Hxy*F(x2D(:),y2D(:));
                    Hk = V'\V';
                    self.b = sem.assembly.assemblybDirac2D(Nk, Np, geofact, self.conn, Hk, self.x2D, self.y2D, source);
                otherwise
                    error('Source type not supported')
            end

            self.L = k^2*self.M - (self.Sx + self.Sy);
            
            [self.Ltilde, self.btilde] = updateBoundaries(self.boundary_type, ...
                self.L,self.b,self.x2D,self.y2D,self.conn,self.acoustic_params,...
                self.mesh_info.bound_detect_f, self.boundary_cond_f);             
        end
        
        function [P,XY] = solveFreq(self, k)
            self.L = k^2*self.M - (self.Sx + self.Sy);
            
            [self.Ltilde, self.btilde] = updateBoundaries(self.boundary_type, ...
                self.L,self.b,self.x2D,self.y2D,self.conn,self.acoustic_params,...
                self.mesh_info.bound_detect_f, self.boundary_cond_f);                    

            [P,XY] = self.solve();
        end

        function [P,XY] = solve(self)            
            self.P = self.Ltilde \ self.btilde;
            P = full(self.P);
            XY = [self.x2D, self.y2D];
        end

        function clearResults(self)
            self.L = 0; self.Ltilde = 0; self.b = 0; self.btilde = 0; self.M = 0; self.Sx = 0; self.Sy = 0;
            self.P = 0; self.x2D = 0; self.y2D = 0;
        end
    end
end

function [Ltilde, btilde] = updateBoundaries(...
    boundary_type,L,b,x2D,y2D,conn,acoustic_params,bound_detect_f,bound_cond_f)

    switch boundary_type
        case models.BoundaryCondition.Velocity
            [Ltilde, btilde] = sem.assembly.assemblyNeumannBounds2D(...
                L,b,x2D,y2D,conn,bound_detect_f,bound_cond_f);
        case models.BoundaryCondition.Pressure
            [Ltilde, btilde] = sem.assembly.assemblyDirichletBounds2D(...
                L,b,x2D,y2D,conn,bound_detect_f,bound_cond_f);
        case models.BoundaryCondition.Impedance
            [Ltilde, btilde] = sem.assembly.assemblyImpedanceBounds2D(...
                L,b,x2D,y2D,conn,acoustic_params,bound_detect_f,bound_cond_f);
        otherwise
            error('boundary type not supported')
    end
end

