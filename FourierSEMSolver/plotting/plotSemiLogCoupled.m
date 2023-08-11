function plotSemiLogCoupled(x1d_1, x1d_2, x1d_ref, p1, p2, p_ref, title_str, legend1, legend2, ax1, ax2)
    t1log = '-or'; t2log = '-ob';   
    t1 = '-ro'; t2 = '-bo'; 
    
    plot(ax1,x1d_1,p1,t1,'LineWidth', 2)
    hold(ax1, 'on')
    plot(ax1,x1d_2,p2,t2,'LineWidth', 2)    
    if isempty(x1d_ref)
        xline(ax1, x1d_1(end), '--')
        legend(ax1, legend1, legend2, 'Interface', 'location', 'southwest')
    else
        plot(ax1,x1d_ref,p_ref, '--')
        xline(ax1, x1d_1(end), '--')
        legend(ax1, legend1, legend2, 'Reference', 'Interface', 'location', 'southwest')
    end    
    
    xlabel(ax1,'x')
    ylabel(ax1,'p')
    title(ax1, title_str)
    ylim(ax1,[-1,1])
    xlim(ax1,[x1d_1(1),x1d_2(end)])    
    set(ax1,'FontSize',20)
    
    plot(ax2,x1d_1,20*log10(abs(p1)),t1log, 'LineWidth', 2) %domain 1
    hold(ax2, 'on')
    plot(ax2,x1d_2,20*log10(abs(p2)),t2log, 'LineWidth', 2) %domain 2
    
    xlabel(ax2,'x')
    ylabel(ax2,'$20\log(|p|)$  [dB]', 'Interpreter', 'latex')
    ylim(ax2,[-80,0])
    xlim(ax2,[x1d_1(1),x1d_2(end)])
    xline(ax2, x1d_1(end), '--')
    legend(ax2, legend1, legend2, 'Interface', 'location', 'southwest')    
    set(ax2,'FontSize',20)
    
    hold(ax1, 'off')
    hold(ax2, 'off')
end