function plotCoupled(state1, state2, p1, p2, scheme_order, label1, label2, ax1)
    t1 = '.r'; t2 = '.b';
        
    x1d_1 = state1.domain.x1d;
    l = state2.domain.xmax - state2.domain.xmin;
    x1d_2 = l + state2.domain.x1d;        

    t = state1.n*state1.domain.dt;
        
    plot(ax1,x1d_1,p1,t1) %domain 1
    hold(ax1, 'on')
    plot(ax1,x1d_2,p2,t2) %domain 2
    
    xlabel(ax1,'x')
    ylabel(ax1,'p')
    title(ax1, sprintf('p_{%s} <--> p_{%s} \n %.4f [s] \n %s order scheme laplacian   dx_1=%.3f   dx_2=%.3f   dt_1=%f   dt_2=%f', .../
        label1, label2, t, scheme_order, state1.domain.dx, state2.domain.dx, state1.domain.dt, state2.domain.dt))
    ylim(ax1,[-1,1])
    %legend(ax1, label1, label2)
    
    hold(ax1, 'off')
end