function [p, etov, XY] = distMeshGeneration(bounding_box, xy, hpoint, hmax)
%
%   hpoint: smallest element around xy will be roughly this size (e.g. 0.01)
%   hmax:   maximum element size (e.g. hmax = 0.25)
%
    xmin = bounding_box(1,1);
    ymin = bounding_box(1,2);

    xmax = bounding_box(2,1);
    ymax = bounding_box(2,2);
    
    c1 = [xmin,ymin];
    c2 = [xmin,ymax];
    c3 = [xmax,ymin];
    c4 = [xmax,ymax];
    
    XY = [c1; c2; c3; c4; xy]; % require source point  

    fh1 = @(p) hpoint + 0.3*dcircle(p,xy(1),xy(2),0);
    fh = @(p) min(fh1(p), hmax);
    fd = @(p) drectangle(p,xmin,xmax,ymin,ymax);
    
    %drectangle = @(p,x1,x2,y1,y2) -min(min(min(-y1+p(:,2),y2-p(:,2)),-x1+p(:,1)),x2-p(:,1));    
    %fh = @(p) 0.05 + 0.3*sqrt(sum(abs((p-xy).^2),2));
    %fh = @(p) hpoint + 0.3*sqrt(sum(abs((p-xy).^1),2));
    %fh = @(p) ones(size(p,1),1);    
    
    [p,etov] = distmesh2d( fd, fh, hpoint, bounding_box, XY);
    
    [p,etov] = fixmesh(p,etov); 
%     % identify unique vertex nodes
%     [p,I,J] = unique(p,'rows');
%     % change numbering of vertex nodes accordingly to sort order
%     etov = reshape(J(etov),size(etov));
%     % remove duplicate entries in p and update EToV accordingly
%     [idx,I,J] = unique(etov);
%     p = p(idx,:);
%     etov = reshape(J,size(etov));
%     
%     N=size(etov,1);         % Number of elements
%     Nv=size(p,1);           % Number of vertex nodes in mesh
%     Nfaces=size(etov,2);    % Number of faces/element
%     VX = p(:,1);            % x-coordinates of vertex nodes
%     VY = p(:,2);            % y-coordinates of vertex nodes
%     
%     etov = meshing.reorder(etov,VX,VY);
    
%     clf
%     figure(10)
%     patch( 'vertices', p, 'faces', etov, 'facecolor', [.9, .9, .9] )
%     title('Mesh')
%     axis tight
%     axis equal
end