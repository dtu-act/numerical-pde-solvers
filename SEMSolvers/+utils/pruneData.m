function [pout, vout, accs_l, accs_r, tout, xout] = pruneData(t, x1d, ...
                                                              p_hat, v_hat, accs_l, accs_r, ...
                                                              xprune, tprune)
    t_pruned = t(1:tprune:end);
    x_pruned = x1d(1:xprune:end);
    
    pout = p_hat(:,1:tprune:end,1:xprune:end);
    vout = v_hat(:,1:tprune:end,1:xprune:end);
    
    accs_l = accs_l(:,1:tprune:end,:);
    accs_r = accs_r(:,1:tprune:end,:);
    
    [xout,tout] = meshgrid(x_pruned, t_pruned);
end