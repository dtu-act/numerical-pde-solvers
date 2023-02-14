function plotConvergenceP(errs, P_orders, Ne_matrix, title_str, path_write)

    if length(P_orders) ~= size(Ne_matrix,1)
        throw("Error: Number of polynomial order tests should match number of rows in the element matrix")
    end

    colors = ['b', 'g', 'r', 'm'];
    
    legendInfo1 = {};
    legendInfo2 = {};
    
    figure(1)
    
    for i=1:size(Ne_matrix,1)
        col = sprintf('%s-o', colors(mod(i-1, length(colors))+1));

        loglog( Ne_matrix(i,:), errs(i,:), col, 'Linewidth',2)
        hold on
        legendInfo1{i} = sprintf('p=%i',P_orders(i));
    end

    for i=1:size(Ne_matrix,1)
        col = sprintf('%s--', colors(mod(i-1, length(colors))+1));
        
        ax = (1 ./ Ne_matrix(i,:).^(P_orders(i)+1));
        b = 2*errs(i,1)/ax(1);
        
        loglog( Ne_matrix(i,:), b*ax, col, 'Linewidth', 2);
        legendInfo2{i} = sprintf('$\\mathcal{O}(\\Delta x^%i)$',P_orders(i)+1);
    end
    hold off
    
    title(title_str)
    legend([legendInfo1, legendInfo2], 'Location','southeast','interpreter','latex')
    xlabel('Number of elements')
    ylabel('L_2 error')
    set(gca,'fontsize',16)
    %ylim([1e-9, 1e-2])
    grid

    if nargin == 5
        saveas(gcf,path_write)
    end
end