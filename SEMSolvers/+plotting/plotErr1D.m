function plotErr1D(p,x,exclude_endpoints, title_str, xlabel_str, ylabel_str)
    extraInputs = {'interpreter','latex','fontsize',14};
    
    if nargin < 4
        xlabel_str = 'x';
        ylabel_str = 'y';
        exclude_endpoints = false;
        title_str = 'Error';
    end
    
    figure()
    if exclude_endpoints
        semilogy(x(2:end-1), p(2:end-1));
    else
        semilogy(x, p);
    end
    
    xlabel(xlabel_str, extraInputs{:})
    ylabel(ylabel_str, extraInputs{:})
    %zlim([0 1.5])
    title(title_str)
    grid
    
    set(gca,'fontsize',18)
end