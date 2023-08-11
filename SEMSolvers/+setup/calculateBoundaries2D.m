function [hsurf, bcmap, M1D] = calculateBoundaries2D(conn, X2D, Y2D, xminmax, yminmax, P, NODETOL)

    xl = xminmax(1);
    xu = xminmax(2);
    yl = yminmax(1);
    yu = yminmax(2);
    
    % Create 1D mass matrix for dealing with BC's
    r1d = spectral.JacobiGL(0,0,P); % local nodes,LGL
    V1D = spectral.Vandermonde1D(P,r1d);
    M1D = inv(V1D*V1D'); % Standard 1D mass for boundary treatment

    % Identify faces that lie on domain boundary
    [nel,nnodepel] = size(conn);
    
    % Counters for # of faces residing on each of the four domain boundary surfaces
    countL = 0;
    countR = 0;
    countT = 0;
    countB = 0;
    
    % Distance functions from arbitrary point to boundaries
    bdyL = @(x0,y0) abs((yu-yl)*x0 - (xl-xl)*y0 + xl*yl - yu*xl)/(sqrt((yu-yl)^2 + (xl-xl)^2));
    bdyR = @(x0,y0) abs((yu-yl)*x0 - (xu-xu)*y0 + xu*yl - yu*xu)/(sqrt((yu-yl)^2 + (xu-xu)^2));
    bdyT = @(x0,y0) abs((yu-yu)*x0 - (xu-xl)*y0 + xu*yu - yu*xl)/(sqrt((yu-yu)^2 + (xu-xl)^2));
    bdyB = @(x0,y0) abs((yl-yl)*x0 - (xu-xl)*y0 + xu*yl - yl*xl)/(sqrt((yl-yl)^2 + (xu-xl)^2));
       
    VNUM = [1 2; 2 3; 3 1]; % Vertex nodes
    NI = P-1; % The number of "interior" edge nodes is P-1
    
    
    for kk = 1:nel % Loop through all elements
        for jj = 1:3 % Loop through all faces within each element
    
            % Left domain boundary
    
            vnodes = [conn(kk,VNUM(jj,1)) conn(kk,VNUM(jj,2))]; % Vertex nodes
            
            dist1 = bdyL(X2D(vnodes(1)), Y2D(vnodes(1))); % Distance from bdy for vertex 1
            dist2 = bdyL(X2D(vnodes(2)), Y2D(vnodes(2))); % Distance from bdy for vertex 2
            
            if dist1 < NODETOL && dist2 < NODETOL % Both vertex nodes of this face are on this boundary
                countL = countL+1;
                % Insert the boundary nodes, in the right order, into the map
                idxs = 4 + (jj-1)*NI; % Find where the "interior" face node start
                idxf = idxs+NI-1; % And where they finish
                inodes = conn(kk,idxs:idxf);% "interior" face nodes
                fnodes = [vnodes(1) inodes vnodes(2)]; % all nodes of the boundary face, with order matching the 1D local mass matrix
                bcLmap(countL,:) = fnodes;
                continue;
            end
    
            % Right domain boundary
            vnodes = [conn(kk,VNUM(jj,1)) conn(kk,VNUM(jj,2))]; % Vertex nodes
    
            dist1 = bdyR(X2D(vnodes(1)), Y2D(vnodes(1))); % Distance from bdy for vertex 1
            dist2 = bdyR(X2D(vnodes(2)), Y2D(vnodes(2))); % Distance from bdy for vertex 2
    
            if dist1 < NODETOL && dist2 < NODETOL % Both vertex nodes of this face are on this boundary
                countR = countR+1;
                % Insert the boundary nodes, in the right order, into the map
                idxs = 4 + (jj-1)*NI; % Find where the "interior" face node start
                idxf = idxs+NI-1; % And where they finish
                inodes = conn(kk,idxs:idxf);% "interior" face nodes
                fnodes = [vnodes(1) inodes vnodes(2)]; % all nodes of the boundary face, with order matching the 1D local mass matrix
                bcRmap(countR,:) = fnodes;
                continue;
            end
            
            % Top domain boundary
            vnodes = [conn(kk,VNUM(jj,1)) conn(kk,VNUM(jj,2))]; % Vertex nodes
    
            dist1 = bdyT(X2D(vnodes(1)), Y2D(vnodes(1))); % Distance from bdy for vertex 1
            dist2 = bdyT(X2D(vnodes(2)), Y2D(vnodes(2))); % Distance from bdy for vertex 2
    
            if dist1 < NODETOL && dist2 < NODETOL % Both vertex nodes of this face are on this boundary
                countT = countT+1;
                % Insert the boundary nodes, in the right order, into the map
                idxs = 4 + (jj-1)*NI; % Find where the "interior" face node start
                idxf = idxs+NI-1; % And where they finish
                inodes = conn(kk,idxs:idxf);% "interior" face nodes
                fnodes = [vnodes(1) inodes vnodes(2)]; % all nodes of the boundary face, with order matching the 1D local mass matrix
                bcTmap(countT,:) = fnodes;
                continue;
            end
    
            
            % Bottom domain boundary
            vnodes = [conn(kk,VNUM(jj,1)) conn(kk,VNUM(jj,2))]; % Vertex nodes
    
            dist1 = bdyB(X2D(vnodes(1)), Y2D(vnodes(1))); % Distance from bdy for vertex 1
            dist2 = bdyB(X2D(vnodes(2)), Y2D(vnodes(2))); % Distance from bdy for vertex 2
    
            if dist1 < NODETOL && dist2 < NODETOL % Both vertex nodes of this face are on this boundary
                countB = countB+1;
                % Insert the boundary nodes, in the right order, into the map
                idxs = 4 + (jj-1)*NI; % Find where the "interior" face node start
                idxf = idxs+NI-1; % And where they finish
                inodes = conn(kk,idxs:idxf);% "interior" face nodes
                fnodes = [vnodes(1) inodes vnodes(2)]; % all nodes of the boundary face, with order matching the 1D local mass matrix
                bcBmap(countB,:) = fnodes;
                continue;
            end
    
        end
    end
    
    countTot = countL+countR+countT+countB; % Total # of faces on the boundary
    
    % Combine boundary maps for each domain boundary (L,R,T,B) to one boundary map
    bcmap = zeros(countTot,P+1);
    count = 0;
    for ii = 1:countL
        count = count+1;
        bcmap(count,:) = bcLmap(ii,:); 
    end
    
    for ii = 1:countR
        count = count+1;
        bcmap(count,:) = bcRmap(ii,:); 
    end
    
    for ii = 1:countT
        count = count+1;
        bcmap(count,:) = bcTmap(ii,:); 
    end
    
    for ii = 1:countB
        count = count+1;
        bcmap(count,:) = bcBmap(ii,:); 
    end
    
    % Calculate length of each of the boundary faces
    hsurf = zeros(1,countTot);
    count = 0;
    for jj = 1:size(bcmap,1)
        count = count+1;
        x1 = X2D(bcmap(jj,1));
        x2 = X2D(bcmap(jj,end));
        y1 = Y2D(bcmap(jj,1));
        y2 = Y2D(bcmap(jj,end));
        
        d = sqrt((x2-x1)^2 + (y2-y1)^2);
        hsurf(count) = d;
    end
end