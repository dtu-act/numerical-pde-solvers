function plotHelmholtzError(P, Pref, X, Y, method_ref, method_val, view_angle)
    assert(length(P) == length(Pref))
    
    if nargin < 7
        view_angle = 90;
    end
    
    if iscell(P)
        N = length(P);
    else
        N = 1;
        P = {P}; Pref = {Pref}; X = {X}; Y = {Y};
    end
    
    for i=1:N
        Pi = P{i};
        Pi_ref = Pref{i};
        
        Xi = X{i};
        Yi = Y{i};

        XYi = [Xi(:), Yi(:)];
        %plotting.plotWBM(P_wbm, X_wbm, Y_wbm)

        figure()
        plotting.plotHelmholtz(abs(Pi(:)), XYi(:,1), XYi(:,2), ...
            sprintf('%s Helmholtz solution, D_%i', method_val, i), [1,3,1], view_angle)
        plotting.plotHelmholtz(abs(Pi_ref(:)), XYi(:,1), XYi(:,2), ...
            sprintf('%s Function solution, D_%i', method_ref, i), [1,3,2], view_angle)
        plotting.plotHelmholtz(abs(abs(Pi_ref(:)) - abs(Pi(:))), XYi(:,1), XYi(:,2), ...
            sprintf('Error, D_%i', i), [1,3,3], view_angle)
    end
end