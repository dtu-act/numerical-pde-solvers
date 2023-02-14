function [conn, x2d, y2d, gidx] = connectivityTable2D(P,X,Y,vx,vy,etov,etoe,etof)
    % output: Tuple{Matrix{Int}, Vector, Vector, Int}
    %
    %  Nv = size(vx,1)       % number of unique (global) vertices
    %

    Nfaces = size(etov,2);     % number of element faces (3 for triangle)
    Nelem = size(etov,1);      % number of total elements (triangles)
    
    gidx = length(vx);         % global index (no assumption on vertices and numbering)
    Mp   = (P + 1)*(P + 2)/2;  % number of points per triangle (?)
    Mpf  = P + 1;              % number of points per triangle element face

    %Nconn = (Nfaces + (Nfaces-1)*(Mpf-2)) + (Mpf-2)
    conn = zeros(Nelem, Mp);

    x2d = vx(sort(unique(etov)));
    y2d = vy(sort(unique(etov)));
    
    for n=1:Nelem
        for i=1:Nfaces

            conn(n,i) = etov(n,i);
            lidx = (Nfaces + (i-1)*(Mpf-2)) + (1:Mpf-2);

            if etoe(n,i) >= n                         % if connecting element is new then we have new nodes                
                conn(n,lidx) = (gidx+1):(gidx+Mpf-2); % insert edge nodes / assign global numbers
                gidx = gidx + Mpf - 2;
                x2d = [x2d; X(lidx,n)];
                y2d = [y2d; Y(lidx,n)];
            else
                % we already assigned global numbers to the edge nodes
                kc = etoe(n,i); % element number of connecting element with edge nodes already set
                ic = etof(n,i); % face number of connecting element       
                
                % reverse ordering in connecting element
                conn(n,lidx) = conn(kc, (Nfaces + (ic - 1)*(Mpf - 2)) + (Mpf-2:-1:1));
            end
        end

        lidx = (Nfaces + Nfaces*(Mpf-2) + 1):Mp;
        
        conn(n,lidx) = (gidx+1):(gidx + Mp - Nfaces - Nfaces*(Mpf - 2));
        gidx = gidx + Mp - Nfaces - Nfaces*(Mpf - 2);

        x2d = [x2d; X(lidx,n)]; % TODO: optimize for speed
        y2d = [y2d; Y(lidx,n)]; % TODO: optimize for speed
    end
end