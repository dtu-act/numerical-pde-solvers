function circles = plotCircle(xy0,r,c)
    hold on
    th = 0:pi/50:2*pi;
    x_circle = r * cos(th) + xy0(1);
    y_circle = r * sin(th) + xy0(2);
    circles = plot(x_circle, y_circle);
    fill(x_circle, y_circle, c)
    hold off    
end