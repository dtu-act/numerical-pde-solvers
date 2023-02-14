function [c,gidx,X2D,Y2D,FaceValues1D] = constructConnectivityTable2D(VX,VY,EToV,P,x,y,varargin)
%
% Construct local-to-global connectivity tables (Algorithm 4.1).
%
% Ordering used is: [Vertices Edgenodes Interiornodes]
%
% By Allan P. Engsig-Karup, apek@dtu.dk.

Nfaces = size(EToV,2);
if nargin>6
    % this one is better when we have slits in the mesh
    EToE = varargin{1};
    EToF = varargin{2};
else
    % just assume that all elements are connected
    [EToE,EToF]= tiConnect2Dquad(EToV); % CHANGED THIS ONE TO QUAD
end
switch Nfaces
    case 3 % triangles
        Np = (P + 1)*(P + 2)/2; % triangles
    case 4 % quadrilaterals
        Np = (P + 1)^2;
    otherwise
        warning('Elements not identified. (triangles/quadrilaterals)')
        wer % force error to break execution.
end
Npf    = P + 1; % number of nodes on an element edge
Nv     = size(VX,1);   % number of global vertices
Nk     = size(EToV,1);
c      = zeros(Nk,Np);
FaceValues1D = zeros(P+1,Nk*Nfaces); % hold global indices of face values for all edges in mesh
X2D    = VX(unique(EToV));
Y2D    = VY(unique(EToV));
gidx   = length(VX);   % no assumption on vertices and numbering
count = 0;
FaceVal = [1 2; 2 3; 3 1];
for k = 1 : Nk
%    disp(sprintf('Element No. %d',k))

    % global vertices
    c(k,1:Nfaces) = EToV(k,:); 

    % global edges
    for i = 1 : Nfaces
        count = count + 1;
%        disp(sprintf('  Adding face No. %d',i))        
        if EToE(k,i)>=k % if connecting element is new then we have new nodes
            % increase global numbering as usual
            lidx = Nfaces+(i-1)*(Npf-2)+(1:Npf-2);
            c(k,lidx) = gidx+(1:Npf-2); % insert edge nodes / assign global numbers
            gidx = gidx + Npf-2;
            tmpx = x(lidx,k); X2D = [X2D tmpx(:)'];
            tmpy = y(lidx,k); Y2D = [Y2D tmpy(:)'];
            
            GidxFace = [c(k,FaceVal(i,1)) c(k,lidx) c(k,FaceVal(i,2))];
            FaceValues1D(:,count) = GidxFace;
        else % we already assigned global numbers to the edge nodes.
            % lookup node numbering
            kconnect = EToE(k,i); % element number of connecting element with edge nodes already set
            iconnect = EToF(k,i);% face number of connecting element
            lidx = Nfaces+(i-1)*(Npf-2)+(1:Npf-2);
            % reverse ordering in connecting element
            c(k,lidx) = c(kconnect,Nfaces+(iconnect-1)*(Npf-2)+(Npf-2:-1:1)); 

%            GidxFace = [c(k,FaceVal(i,1)) c(k,lidx(end:-1:1)) c(k,FaceVal(i,2))];
            GidxFace = [c(k,FaceVal(i,1)) c(k,lidx) c(k,FaceVal(i,2))];
            FaceValues1D(:,count) = GidxFace;
        end
    end
%    disp(sprintf('  Adding interior nodes'))

    % global interior modes
    lidx = Nfaces+Nfaces*(Npf-2)+1:Np; % local indices for interior nodes
    c(k,lidx) = gidx + (1 : Np - Nfaces - Nfaces*(Npf-2) ); % add new interior nodes
    tmpx = x(lidx,k); X2D = [X2D tmpx(:)'];
    tmpy = y(lidx,k); Y2D = [Y2D tmpy(:)'];    
    gidx = gidx + Np - Nfaces - Nfaces*(Npf-2);
end
%disp(sprintf('Total No. Elements is %d.',Nk))
return




%% Show face nodes
figure
triplot(EToV,VX,VY,'k')
hold on
%plot(X2D,Y2D,'k.')
FaceValues = zeros(N+1,Nk2D*Nfaces);
Npf = N+1;
count = 0;
for k = 1 : Nk2D
triplot(EToV,VX,VY,'k')
hold on
    count = count + 1;
    lidx = [1 Nfaces+(2:Npf-1) 2];
    GidxF1 = c2D(k,lidx);
    plot(X2D(GidxF1),Y2D(GidxF1),'rd')
    FaceValues(:,count) = GidxF1;
    count = count + 1;
    lidx = [2 Nfaces+(Npf-2)+(2:Npf-1) 3];
    GidxF2 = c2D(k,lidx);
    plot(X2D(GidxF2),Y2D(GidxF2),'go')
    FaceValues(:,count) = GidxF2;
    count = count + 1;
    lidx = [3 Nfaces+2*(Npf-2)+(2:Npf-1) 1];
    GidxF3 = c2D(k,lidx);
    plot(X2D(GidxF3),Y2D(GidxF3),'mx')
    FaceValues(:,count) = GidxF3;
    hold off
    drawnow
%    pause
end


    