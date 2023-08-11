function plotConvergenceMesh(errs, xvalues, Proi, XY, xy0, xlabelstr, legendstr, titlestr, ...
    line_style, angle_plot)
    
    subplot(1,2,2)
    loglog( xvalues, errs, line_style, 'Linewidth', 2)
    
    hold on
    
    Porders_ref = [1,2,4,6];
    line_styles_ref = ["--b", "--k", "--r", "--m"];
    
    legendInfo2 = {};
    legendInfo1 = {legendstr};

    for i=1:length(Porders_ref)
        ax = (1 ./ xvalues.^(Porders_ref(i)+1));
        b = errs(1)/ax(1);
        hold on
        loglog( xvalues, b*ax, line_styles_ref(i), 'Linewidth', 2);
        legendInfo2{i} = sprintf('$\\mathcal{O}(\\Delta x^%i)$',Porders_ref(i)+1);
    end
    xlabel(xlabelstr)
    ylabel('L_{2} error')
    title(titlestr)    
    legend([legendInfo1, legendInfo2], 'Location','southeast','interpreter','latex')
    set(gca,'fontsize',18)
    grid

    plotting.plotHelmholtz(abs(Proi(:)), XY(:,1), XY(:,2), legendstr, [1,2,1], angle_plot)
    plotting.plotCircle(xy0,0.3,'w');
    hold off
end