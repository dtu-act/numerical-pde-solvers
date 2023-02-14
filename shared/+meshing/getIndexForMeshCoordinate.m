function [i] = getIndexForMeshCoordinate(XY,x0,y0)
    irows = find(abs(XY(:,1)-x0) == min(abs(XY(:,1)-x0)));
    [~,j] = min(abs(XY(irows,2) - y0));

    i = irows(j);
end