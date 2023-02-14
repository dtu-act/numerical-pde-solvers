function plotHelmholtzDomains(P, X, Y)
    N = length(P);
    
    for i=1:N
        Pi = P{i};
        
        Xi = X{i};
        Yi = Y{i};

        XYi = [Xi(:), Yi(:)];
        %plotting.plotWBM(P_wbm, X_wbm, Y_wbm)

        figure(1)
        plotting.plotHelmholtz(abs(Pi(:)), XYi(:,1), XYi(:,2), sprintf('WBM Helmholtz solution, D_%i', i), [1,N,i])

        figure(2)
        if i == 1
            % HACK! 
            Xi = flip(Xi,1);
            Yi = flip(Yi,1);
        end
        mesh( Xi(1,:), Yi(:,1), abs(Pi), 'linewidth', 1);

        hold on
        shading flat

        view(0,90)
        set(gcf, 'WindowState', 'maximized');
        colorbar
    end
end