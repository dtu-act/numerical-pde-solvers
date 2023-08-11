function writeGif(handle, path, init)    
    % Capture the plot as an image 
    frame = getframe(handle); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if init
        imwrite(imind,cm,path,'gif', 'Loopcount',inf); 
    else
        imwrite(imind,cm,path,'gif','DelayTime',0.005,'WriteMode','append');
    end 
end