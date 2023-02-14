function [M2D,Sx,Sy] = assembleWaveEquationCoupled2D(Nk,Np,conn,x,y,r,s,P)
    V2D = spectral.Vandermonde2D(P,r,s);
    [Dr2D,Ds2D] = spectral.Dmatrices2D(P, r, s, V2D);
    gfs = spectral.GeometricFactors2D(x,y,Dr2D,Ds2D);
    M  = inv(V2D*V2D');

    %% GLOBAL ASSEMBLY
    [M2D] = globalAssemblySEMMass(Nk,Np,conn,M,gfs.J);
    [Sx]  = globalAssemblySEMStiffness(Nk,Np,conn,M,gfs.J,gfs.rx,gfs.sx,Dr2D,Ds2D);
    [Sy]  = globalAssemblySEMStiffness(Nk,Np,conn,M,gfs.J,gfs.ry,gfs.sy,Dr2D,Ds2D);
end

function [A] = globalAssemblySEMMass(Nk,Np,conn,M,J2D)
    Ai = ones(Np*Np*Nk,1);
    Aj = ones(Np*Np*Nk,1);
    As = zeros(Np*Np*Nk,1);
    count = 0;

    for k = 1 : Nk
        Jk = J2D(:,k);
        Ak = diag(Jk)*M;
        idx = conn(k,:);
        for j = 1 : Np
            Ai(count + (1:Np)) = idx;
            Aj(count + (1:Np)) = idx(j);
            As(count + (1:Np)) = Ak(:,j);
            count = count + Np;
        end
    end
    
    A = sparse(Ai,Aj,As);
end

function [A] = globalAssemblySEMStiffness(Nk,Np,conn,M,J2D,rx2D,sx2D,Dr2D,Ds2D)
    % {rx | sx | ry | sy | J | xr | xs | yr | ys} 

    Ai = ones(Np*Np*Nk,1);
    Aj = ones(Np*Np*Nk,1);
    As = zeros(Np*Np*Nk,1);
    count = 0;

    for k = 1:Nk
        DXk = diag(rx2D(:,k))*Dr2D+diag(sx2D(:,k))*Ds2D;  
        % Local assembly term by term
        Jk = J2D(:,k);
        Ak = diag(Jk)*M*DXk;% Mass matrix, no b
        idx = conn(k,:);
        for j = 1 : Np
            Ai(count + (1:Np)) = idx;
            Aj(count + (1:Np)) = idx(j);
            As(count + (1:Np)) = Ak(:,j);
            count = count + Np;
        end
    end
    
    A = sparse(Ai,Aj,As);
end