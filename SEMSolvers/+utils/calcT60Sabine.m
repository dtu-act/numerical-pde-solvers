function T_60 = calcT60Sabine(bbox,alpha)
    lx = bbox(2,1) + bbox(1,1); ly = bbox(2,2) + bbox(1,2);
    V = lx*ly; % volume
    S = lx*ly; % surface
    T_60=0.16*V/(alpha*S); % sabine
end