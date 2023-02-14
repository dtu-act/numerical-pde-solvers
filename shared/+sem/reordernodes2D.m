function [x_reordered,y_reordered,r_reordered,s_reordered] = reordernodes2D(P,x,y,r,s)
%
%    Reorder as: vertex nodes first, then edge nodes and last all interior nodes
%
%    See p. 81 "The Spectral/hp-Finite Element Method for PDE" by A. Engsig-Karup
%

    % Compute index maps for node positions
    tol = min(0.2/P^2,1e-5);
    Mp = length(s);
    
    % Faces
    fid1 = find( abs(s+1) < tol)';
    fid2 = find( abs(r+s) < tol)';
    fid3 = find( abs(r+1) < tol)';
    
    % Interior
    fint = sort(setdiff(1:Mp,[fid1 fid2 fid3]));
    
    % Reorder local nodes to match the chosen ordering in terms of vertex nodes,
    % edge nodes and interior nodes.
    Mpf = P+1;

    localReorder = [1 Mpf Mp fid1(2:Mpf-1) fid2(2:Mpf-1) fid3(Mpf-1:-1:2) fint];
    
    if length(localReorder) ~= length(r(:))
        throw("tol is not correct");
    end

    x_reordered = x(localReorder,:);
    y_reordered = y(localReorder,:);
    r_reordered = r(localReorder);
    s_reordered = s(localReorder);
end