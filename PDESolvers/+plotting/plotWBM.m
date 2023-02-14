function plotWBM(P, rx, ry)
    disp("Plotting....\n")

    xs = rx(1,:);
    ys = ry(:,1)';

    xs1d = rx(:);
    ys1d = ry(:);

    tri = delaunay(xs1d,ys1d);

    figure(1)
    contourf(xs, ys, abs(P));
    colorbar;
    xlabel('x')
    ylabel('y')
    title('WBM Helmholtz solution')
    figure(2)
    trisurf(tri,xs1d,ys1d,abs(P));
    colorbar;
    xlabel('x')
    ylabel('y')
    title('WBM Helmholtz solution')
%     figure(3)
%     trisurf(tri,xs1d,ys1d,P_err);
%     title('Err')
end