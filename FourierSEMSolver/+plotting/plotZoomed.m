function plotZoomed(n, x1d, p1, iface_order, nindx, title_str, legends_str)
    colororder(['b'; 'r'])
    
    if mod(n, 5) == 0
        subplot(2,1,1)        
        plot(x1d,p1,'-')        
        hold on
        %plot(x1d(2:iface_order/2+1,:),p1(2:iface_order/2+1,:),"or")
        hold off
        xlabel('x')
        ylabel('p')
        ylim([-1,1])
        xlim([min(x1d(1,:)),max(x1d(end,:))])       
        title(title_str)
        if nargin == 9
            legend(legends_str)
        end

        subplot(2,1,2)
        plot(x1d(1:nindx,:),p1(1:nindx,:),'-o')
        hold on
        plot(x1d(2:iface_order/2+1,:),p1(2:iface_order/2+1,:),"or")
        hold off
        xlabel('x')
        ylabel('p')
        ylim([-1,1])
        xlim([min(x1d(1,:)),max(x1d(nindx,:))])
        title('ZOOM left')
        if nargin == 9
            legend(legends_str)
        end

        drawnow
    end
end