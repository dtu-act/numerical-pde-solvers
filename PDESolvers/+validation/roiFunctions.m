function fd = roiFunctions(roi_type,bthickness,slice_source,bbox,xy0,r0)
%% ROI types
% 'fd_leftb','fd_rightb','fd_lowerb','fd_upperb', 'fd_boundaries'
% 'fd_xaxial','fd_yaxial','fd_axial',
% 'fd_Q1','fd_Q2','fd_Q3','fd_Q4','fd_Q1Q2','fd_Q3Q4','fd_Q1Q2Q3Q4'
% 'fd_bsource', 'fd_bsourcenoaxial'
% 'fd_full', 'fd_fullnobounds'

    assert(bbox(1,1) == 0);
    assert(bbox(1,2) == 0);
    
    lx = bbox(2,1);
    ly = bbox(2,2);

    %% BOUNDARIES
    xmin_leftb = 0.0; xmax_leftb = bthickness; ymin_leftb = 0.0; ymax_leftb = ly;
    xmin_rightb = ly-bthickness; xmax_rightb = lx; ymin_rightb = 0.0; ymax_rightb = ly;
    xmin_lowerb = 0.0; xmax_lowerb = lx; ymin_lowerb = 0.0; ymax_lowerb = bthickness;
    xmin_upperb = 0.0; xmax_upperb = lx; ymin_upperb = ly-bthickness; ymax_upperb = ly;    
    
    fd_leftb = @(p) min(-drectangle(p,xmin_leftb,xmax_leftb,ymin_leftb,ymax_leftb), dcircle(p,xy0(1),xy0(2),r0));
    fd_rightb = @(p) min(-drectangle(p,xmin_rightb,xmax_rightb,ymin_rightb,ymax_rightb), dcircle(p,xy0(1),xy0(2),r0));
    fd_lowerb = @(p) min(-drectangle(p,xmin_lowerb,xmax_lowerb,ymin_lowerb,ymax_lowerb), dcircle(p,xy0(1),xy0(2),r0));
    fd_upperb = @(p) min(-drectangle(p,xmin_upperb,xmax_upperb,ymin_upperb,ymax_upperb), dcircle(p,xy0(1),xy0(2),r0));
    fd_boundaries = @(p) max(max(max(fd_leftb(p), fd_rightb(p)), fd_lowerb(p)), fd_upperb(p));
    
    %% AXIAL 
    xmin_xaxial = bthickness; xmax_xaxial = lx-bthickness; ymin_xaxial = xy0(1)-bthickness/2; ymax_xaxial = xy0(1)+bthickness/2;
    xmin_yaxial = xy0(1)-bthickness/2; xmax_yaxial = xy0(1)+bthickness/2; ymin_yaxial = bthickness; ymax_yaxial = ly-bthickness;

    fd_xaxial = @(p) min(-drectangle(p,xmin_xaxial,xmax_xaxial,ymin_xaxial,ymax_xaxial), dcircle(p,xy0(1),xy0(2),r0));
    fd_yaxial = @(p) min(-drectangle(p,xmin_yaxial,xmax_yaxial,ymin_yaxial,ymax_yaxial), dcircle(p,xy0(1),xy0(2),r0));
    fd_axial = @(p) max(fd_xaxial(p),fd_yaxial(p));

    %% QUADRANTS   
    xmin_Q1 = xy0(1)+bthickness; xmax_Q1 = lx-bthickness; ymin_Q1 = xy0(2)+bthickness; ymax_Q1 = ly-bthickness;
    xmin_Q2 = bthickness; xmax_Q2 = xy0(1)-bthickness; ymin_Q2 = xy0(2)+bthickness; ymax_Q2 = ly-bthickness;
    xmin_Q3 = bthickness; xmax_Q3 = xy0(1)-bthickness; ymin_Q3 = bthickness; ymax_Q3 = xy0(2)-bthickness;
    xmin_Q4 = xy0(1)+bthickness; xmax_Q4 = lx-bthickness; ymin_Q4 = bthickness; ymax_Q4 = xy0(2)-bthickness;

    fd_Q1 = @(p) min(-drectangle(p,xmin_Q1,xmax_Q1,ymin_Q1,ymax_Q1), dcircle(p,xy0(1),xy0(2),r0));
    fd_Q2 = @(p) min(-drectangle(p,xmin_Q2,xmax_Q2,ymin_Q2,ymax_Q2), dcircle(p,xy0(1),xy0(2),r0));
    fd_Q3 = @(p) min(-drectangle(p,xmin_Q3,xmax_Q3,ymin_Q3,ymax_Q3), dcircle(p,xy0(1),xy0(2),r0));
    fd_Q4 = @(p) min(-drectangle(p,xmin_Q4,xmax_Q4,ymin_Q4,ymax_Q4), dcircle(p,xy0(1),xy0(2),r0));

    fd_Q1Q2 = @(p) max(fd_Q1(p), fd_Q2(p));
    fd_Q3Q4 = @(p) max(fd_Q3(p), fd_Q4(p));
    fd_Q1Q2Q3Q4 = @(p) max(fd_Q1Q2(p), fd_Q3Q4(p));

    %% CIRCULAR
    fd_bsource = @(p) min(-dcircle(p,xy0(1),xy0(2),r0+slice_source), dcircle(p,xy0(1),xy0(2),r0));
    fd_bsourcenoaxial = @(p) min(min(-dcircle(p,xy0(1),xy0(2),r0+slice_source),...
                                   dcircle(p,xy0(1),xy0(2),r0)), -fd_axial(p) );

    %% FULL DOMAINS
    fd_full = @(p) min(-drectangle(p,0,lx,0,ly), dcircle(p,xy0(1),xy0(2),r0));
    fd_fullnobounds = @(p) min(-drectangle(p,bthickness,lx-bthickness,bthickness,ly-bthickness), dcircle(p,xy0(1),xy0(2),r0));
    
    switch roi_type
        case 'fd_leftb'
            fd = fd_leftb;
        case 'fd_rightb'
            fd = fd_rightb;
        case 'fd_lowerb'
            fd = fd_lowerb;
        case 'fd_upperb'
            fd = fd_upperb;
        case 'fd_boundaries'
            fd = fd_boundaries;
            
        case 'fd_xaxial'
            fd = fd_xaxial;
        case 'fd_yaxial'
            fd = fd_yaxial;
        case 'fd_axial'
            fd = fd_axial;
            
        case 'fd_Q1'
            fd = fd_Q1;
        case 'fd_Q2'
            fd = fd_Q2;
        case 'fd_Q3'
            fd = fd_Q3;
        case 'fd_Q4'
            fd = fd_Q4;

        case 'fd_Q1Q2'
            fd = fd_Q1Q2;
        case 'fd_Q3Q4'
            fd = fd_Q3Q4;
        case 'fd_Q1Q2Q3Q4'
            fd = fd_Q1Q2Q3Q4;
            
        case 'fd_bsource'
            fd = fd_bsource;
        case 'fd_bsourcenoaxial'
            fd = fd_bsourcenoaxial;
           
        case 'fd_fullnobounds'
            fd = fd_fullnobounds;
        case 'fd_full'
            fd = fd_full;
           
    end

end