function C = calcNormConstantCoord(P1,P2,XY,coords)
    C = 0;
    for k=1:size(coords,1)            
        x = coords(k,1);
        y = coords(k,2);            

        irows = find(abs(XY(:,1)-x) == min(abs(XY(:,1)-x)));
        [~,j] = min(abs(XY(irows,2) - y));

        i = irows(j);
        C = C + P2(i)/P1(i);
    end
    C = C/size(coords,1);
end