function plotTemporalInterpolations(ax1, x1d, state_full, state_exact, interpolations)

    hold(ax1, 'on')
    plot(ax1, x1d, state_full.p_current, 'o-r');    
    plot(ax1, x1d, state_full.p_prev, 'o-b');
    plot(ax1, x1d, state_exact.p_current, '--k');
    % nhistory = 5, hence we get 9 interpolations in total
    %plot(ax1, x1d, interpolations(:,1));
    plot(ax1, x1d, interpolations(:,2));
    %plot(ax1, x1d, interpolations(:,3));
    %plot(ax1, x1d, interpolations(:,4));
    %plot(ax1, x1d, interpolations(:,5));
    %plot(ax1, x1d, interpolations(:,6));
    %plot(ax1, x1d, interpolations(:,7));
    %plot(ax1, x1d, interpolations(:,8));
    %plot(ax1, x1d, interpolations(:,9));
    xlabel(ax1, 'x')
    ylabel(ax1, 'p')
    title(ax1, 'Temporal spline interpolations')
    legend(ax1, 'p^{n+1}', 'p^{n}', 'p_{exact}^{n+1/2}', 'p_{interpolated}^{n+1/2}')
    hold(ax1, 'off')
end