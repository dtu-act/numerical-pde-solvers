function plotConvergences(errs_all, xvalues_all, Porders_plot, xlabelstr, legendInfo1, ...
    line_styles, title_str)
    
    line_styles_ref = ["--b", "--k", "--r", "--m"];
    
    legendInfo2 = {};
    
    for i=1:size(errs_all,1)        
        loglog( xvalues_all{i}, errs_all{i}, line_styles(i), 'Linewidth',2)
        hold on
    end

    for i=1:length(Porders_plot)
        ax = (1 ./ xvalues_all{1}.^(Porders_plot(i)+1));
        b = 2*errs_all{1}(1)/ax(1);
        
        loglog( xvalues_all{1}, b*ax, line_styles_ref(i), 'Linewidth', 2);
        legendInfo2{i} = sprintf('$\\mathcal{O}(\\Delta x^%i)$',Porders_plot(i)+1);
    end
    hold off
    
    title(title_str)
    legend([legendInfo1, legendInfo2], 'Location','southwest','interpreter','latex')
    xlabel(xlabelstr)
    ylabel('L_2 error')
    set(gca,'fontsize',16)
    %ylim([1e-9, 1e-2])
    grid
    set(gca,'xticklabel',arrayfun(@(x) num2str(x),get(gca,'xtick'),'un',0))
    
    set(gcf,'color','w','position',[200 200 600 700]);
end