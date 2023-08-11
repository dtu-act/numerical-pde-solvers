function visualizeMesh(XY, etov, conn)
    x2D = XY(:,1);
    y2D = XY(:,2);
    
    if size(conn,2) > 3
        conn_to_etov = delaunay(x2D,y2D);
        
        figure()
        subplot(1,2,1)        
        visualizeOnTop(XY, etov, etov);
        title('Grid points')
        subplot(1,2,2)        
        visualizeOnTop(XY, conn_to_etov, conn);
        title('All nodes')
    else
        figure()      
        visualizeOnTop(XY, etov, conn);
    end
    
    
end

function visualizeOnTop(XY, etov, conn)
    x2D = XY(:,1);
    y2D = XY(:,2);
    
    nodal_indexes = [];    
    simpplot(XY, etov)
    hold on
    
    for i = 1:size(conn,1)
        for j=conn(i,:)
            if ~ismember(j, nodal_indexes)
                nodal_indexes = [nodal_indexes j];
                scatter(x2D(j),y2D(j), 'k', 'filled')
                b = num2str(j); ccc = cellstr(b);
                dx = 0.02; dy = 0.02; % displacement so the text does not overlay the data points
                text(x2D(j)+dx, y2D(j)+dy, ccc);
            end
        end
    end
    
    hold off
end