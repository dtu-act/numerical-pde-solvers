function geometryPlot(domains, fig_i)
    N_domains = length(domains);
    
    if N_domains == 1 && ~iscell(domains)
        domains = {domains};
    end
    
    colors = ['b', 'g', 'r', 'm'];
    
    if nargin == 1
        figure()
    else
        figure(fig_i)        
    end
    
    xymin = 0;
    xymax = 0;
    
    for i=1:N_domains
        domain = domains{i};
        col = sprintf('%s', colors(mod(i-1, length(colors))+1));
        
        gcoords = domain.gcoord;
        boundary_nodes = domain.nodes_bound;        

        xymin_i = min(gcoords);
        xymax_i = max(gcoords);
        
        centerpoint_XY =[(xymax_i(1)+xymin_i(1))/2, (xymax_i(2)+xymin_i(2))/2];
        text(centerpoint_XY(1),centerpoint_XY(2), sprintf('D_%i', i),'Color','k','FontSize',18)
        
        xymin = min(xymin_i, xymin);
        xymax = max(xymax_i, xymax);
            
        %% Drawing geometry; First: outerboundary, Second; interface boundary
        for k = 1: length(boundary_nodes)
            x1=gcoords(boundary_nodes(k,1),1);
            y1=gcoords(boundary_nodes(k,1),2);
            x2=gcoords(boundary_nodes(k,2),1);
            y2=gcoords(boundary_nodes(k,2),2);

            P_x = [x1,x2];  P_y = [y1,y2];
            line(P_x,P_y,'Color','k','LineWidth',2)

            % Outer Line numbering
            centerpoint_XY =[(x2+x1)/2, (y2+y1)/2];
            pxy = paddingXY(gcoords, centerpoint_XY);
            text(centerpoint_XY(1)+pxy(1),centerpoint_XY(2)+pxy(2), sprintf('%i_%i', k, i),'Color','r','FontSize',14)
            % Node marker drawing
            hold on
            plot(P_x,P_y,'ro','MarkerEdgeColor','b','MarkerSize',8);
        end
        
        bound_id = k;
        
        interface_nodes = domain.nodes_interface;
        N_interfaceNode = size(interface_nodes);

        if ~isempty(N_interfaceNode)      
            for j = 1: N_interfaceNode(1)            
                x1=gcoords(interface_nodes(j,1),1);
                y1=gcoords(interface_nodes(j,1),2);
                x2=gcoords(interface_nodes(j,2),1);
                y2=gcoords(interface_nodes(j,2),2);

                P_x = [ x1,x2];  P_y = [ y1,y2];
                line(P_x,P_y,'Color','g','LineWidth',2,'LineStyle','-.')
                hold on
                plot(P_x,P_y,'ro','MarkerEdgeColor','b','MarkerSize',8);
                centerpoint_XY =[(x2+x1)/2, (y2+y1)/2 ];

                pxy = paddingXY(gcoords, centerpoint_XY);
                text(centerpoint_XY(1)+pxy(1),centerpoint_XY(2)+pxy(1), sprintf('%i_%i', bound_id+j, i),'Color','r','FontSize',14)
            end
        end
        
        grid on
        %% Node number labeling
        for k = 1:length(gcoords)
            pxy = paddingXY(gcoords, [gcoords(k,1),gcoords(k,2)]);
            text(gcoords(k,1)+pxy(1),gcoords(k,2)+pxy(2), sprintf('%i_%i', k, i), 'Color','b','FontSize',14)
        end
    end
    
    set(gcf, 'WindowState', 'maximized');
    axis( [xymin(1)-0.11 xymax(1)+0.11 xymin(2)-0.11 xymax(2)+0.11 ])
    axis equal
    hold off
end

function padding_xy = paddingXY(gcoords,gcoord)
    padding_xy = [0,0];
    if min(gcoords(:,1)) == gcoord(1)
        padding_xy(1) = 0.2;
    elseif max(gcoords(:,1)) == gcoord(1)
        padding_xy(1) = -0.2;
    end
    
    if min(gcoords(:,2)) == gcoord(2)
        padding_xy(2) = 0.2;
    elseif max(gcoords(:,2)) == gcoord(2)
        padding_xy(2) = -0.2;
    end
end
