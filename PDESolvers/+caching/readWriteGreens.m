function [Pref] = readWriteGreens(path,bbox,boundary_type,XY,xy0,aparams,Nw_greens)
    if isfile(path)
        fprintf('[cache] Loading Greens function from disk ... \n')
        data = load(path, 'P');
        Pref = data.P;
    else
        fprintf('[cache] Calculating Greens function ... \n')

        switch boundary_type
            case models.BoundaryCondition.Pressure
                P = greens.greensDirichletUnstructuredSpatial2D(Nw_greens,Nw_greens,XY,xy0,aparams.k);
            case models.BoundaryCondition.Velocity
                P = greens.greensNeumannUnstructuredSpatial2D(Nw_greens,Nw_greens,XY,xy0,aparams.k);
            case models.BoundaryCondition.Impedance
                % alpha = 1 - ((Z - rho*c)/(Z - rho*c))^2
                alpha = utils.ZtoAlpha(aparams.Z,aparams.c,aparams.rho);
                T_60 = utils.calcT60Sabine(bbox,alpha);
                tau = T_60/13.8;
                P = greens.greensNeumannUnstructuredSpatial2D(Nw_greens,Nw_greens,XY,xy0,...
                    aparams.k,tau,aparams.c);
            otherwise
                error('boundary type not supported')
        end

        fprintf('[cache] Writing Greens function to disk: %s... \n', path)            
        save(path, 'P');
        Pref = P;
    end
end