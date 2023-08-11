function plotHelmholtz(P, rx, ry, title_str, subplot_ijk, view_angle)
    disp("Plotting....\n")

    xs1d = rx(:);
    ys1d = ry(:);

    tri = delaunay(xs1d,ys1d);
        
    extraInputs = {'interpreter','latex','fontsize',14};
    
    if nargin < 6
        view_angle_x = 0;
        view_angle_y = 90;        
    else
        view_angle_x = view_angle;
        view_angle_y = view_angle;
    end
    
    if nargin == 3
        title_str = "Helmholz solution";
    end
    
    if nargin >= 5
        subplot(subplot_ijk(1),subplot_ijk(2),subplot_ijk(3))    
    else
        figure()
    end
    trisurf(tri,xs1d,ys1d, P);
    colorbar;
    xlabel('x', extraInputs{:})
    ylabel('y', extraInputs{:})
    %zlim([0 0.4])
    title(title_str)
    
    set(gca,'fontsize',18)
    view(view_angle_x,view_angle_y)
    set(gcf,'color','w','position',[200 200 1200 700]);
    h = get(gca,'DataAspectRatio');
    if h(3)==1
        set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))])
    else
        set(gca,'DataAspectRatio',[1 1 h(3)])
    end
    set(gcf, 'WindowState', 'maximized');
    
%     hold on
%     xmin = 0.2; xmax = 1.0; ymin = 0.2; ymax = 1.0;
%     h = fill3([xmin, xmax,xmax,xmin],[ymin,ymin,ymax,ymax],[1,1,1,1], 'g');
%     h.FaceAlpha=0.3;
%     hold off
end