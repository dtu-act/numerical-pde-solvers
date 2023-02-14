function xy0 = calculateSourcePosition(NeX,NeY,bbox,xy_rel_pos)
    assert(xy_rel_pos(1) <= 1 && xy_rel_pos(2) <= 1);
    
    Lx = (bbox(1,1)+bbox(2,1));
    Ly = (bbox(1,2)+bbox(2,2));
    dx = Lx/(NeX);
    dy = Ly/(NeY);
    
    x0 = round(NeX*xy_rel_pos(1))*dx;
    y0 = round(NeY*xy_rel_pos(2))*dy;
    xy0 = [x0,y0];
end