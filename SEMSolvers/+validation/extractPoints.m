function indxs = extractPoints(XY,xys)
    N = size(xys,1);
    indxs = zeros(1,N);
    for i=1:N
        indxs(i) = meshing.getIndexForMeshCoordinate(XY,xys(i,1),xys(i,2));        
    end
    indxs = unique(indxs);
end