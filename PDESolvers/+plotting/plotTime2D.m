function plotTime2D(F_time,Fref_time,T,XY,rec_point,plot_frac,file_path,title_str)

    N = size(F_time,1);
    dt = T/N;
    
    if nargin <= 5
        plot_frac = 1;        
    end    
    if nargin <= 6
        file_path = {};
    end
    if nargin <= 7
        title_str = 'Time plot 2D';
    end
    
    disp("Plotting....\n")

    X = XY(:,1);
    Y = XY(:,2);
    
    indx_receiver = meshing.getIndexForMeshCoordinate(XY,rec_point(1),rec_point(2));
    x0 = XY(indx_receiver,1);
    y0 = XY(indx_receiver,2);
    
    tri = delaunay(X,Y);
    
    zmin = min(F_time, [], 'all');
    zmax = max(F_time, [], 'all');
    
    ymin = min(F_time(:,indx_receiver));
    ymax = max(F_time(:,indx_receiver));
    
    ts = []; Fs = []; Fs_ref = [];
    
    h = figure();
    sg = sgtitle(title_str);
    sg.FontSize = 18;
    sg.FontWeight = 'bold';
    
    set(gcf, 'Position',  [100, 100, 1200, 900])
    for i = 1:N
        
        if mod(i,plot_frac) == 0
            t = i*dt;
            
            subplot(2,1,1)
            trisurf(tri,X,Y,F_time(i,:));
            xlabel('x [m]')
            ylabel('y [m]')
            zlabel('p')
            title(sprintf('SEM, t = %.4f [sec]',t))
            zlim([zmin,zmax])
            set(gca,'fontsize',16)
            
            subplot(2,1,2)
            ts = [ts t];
            Fs = [Fs F_time(i,indx_receiver)];
            Fs_ref = [Fs_ref Fref_time(i,indx_receiver)];            
            plot(ts, Fs, '-b');
            hold on
            plot(ts, Fs_ref, '--r');
            hold off            
            title(sprintf('IR at (r_x,r_y) = (%0.2f,%0.2f)', x0,y0))
            xlim([1*dt,N*dt]);
            %ylim([ymin,ymax]);
            xlabel('t [sec]')
            ylabel('Amplitude')
            legend('SEM', 'Analytical')
            set(gca,'fontsize',16)
            
            drawnow
            
            if ~isempty(file_path)
                writeGif(i,h,file_path)
            end
            pause(0.001);
        end
    end
end

% https://se.mathworks.com/help/matlab/ref/imwrite.html
function writeGif(n,h,filename)
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if n == 1 
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.05);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
    end 
end