function [ps_all, vs_all, laccs_all, raccs_all] = rk4_coupled_weq(dxV, dxP, ps0, vs0, t0, tlast, timesteps, dtAccs)
    if length(ps0) ~= length(vs0)
        error("The coupled system doesn't have same spatial size")
    end
    
    if nargin == 7
        dtAccs = @(accs, p) 0;
    end

    ps_all = zeros(timesteps, length(ps0));
    ps_all(1,:) = ps0;
    vs_all = zeros(timesteps, length(ps0));
    vs_all(1,:) = vs0;    

    ps = zeros(size(ps0,1),1);
    ps(:) = ps0;
    vs = zeros(size(ps0,1),1);
    vs(:) = vs0;
        
    num_accs = 4;
    laccs = zeros(num_accs,1);
    raccs = zeros(num_accs,1);
    laccs_all = zeros(timesteps,num_accs);
    raccs_all = zeros(timesteps,num_accs);

    h = (tlast-t0)/timesteps;
    
    t = 0.0;
    
    for i=1:timesteps-1
        % calculate the pressure and the velocity for the next timestep

        K1v = h * dxV(t,          ps);      % calculation of velocity is dependent on pressure
        K1p = h * dxP(t,          vs, ps, laccs, raccs);      % calculation of pressure is dependent on velocity
        K1laccs = h * dtAccs(laccs, ps(1));   % update accumulators
        K1raccs = h * dtAccs(raccs, ps(end)); % update accumulators

        K2v = h * dxV(t + h/2,  ps + K1p/2);
        K2p = h * dxP(t + h/2,  vs + K1v/2, ps + K1p/2, laccs + K1laccs/2, raccs + K1raccs/2);
        K2laccs = h * dtAccs(laccs + K1laccs, ps(1) + K1p(1)/2);
        K2raccs = h * dtAccs(raccs + K1raccs, ps(end) + K1p(end)/2);

        K3v = h * dxV(t + h/2,  ps + K2p/2);
        K3p = h * dxP(t + h/2,  vs + K2v/2, ps + K2p/2, laccs + K2laccs/2, raccs + K2raccs/2);
        K3laccs = h * dtAccs(laccs + K2laccs, ps(1) + K1p(1)/2);
        K3raccs = h * dtAccs(raccs + K2raccs, ps(end) + K1p(end)/2);

        K4v = h * dxV(t + h,      ps + K3p);
        K4p = h * dxP(t + h,      vs + K3v, ps + K3p, laccs + K3laccs, raccs + K3raccs);
        K4lacc = h * dtAccs(laccs + K3laccs, ps(1) + K1p(1)/2);
        K4racc = h * dtAccs(raccs + K3raccs, ps(end) + K1p(end)/2);

        vs = vs + 1/6 * (K1v + 2*(K2v + K3v) + K4v);
        ps = ps + 1/6 * (K1p + 2*(K2p + K3p) + K4p);
        laccs = laccs + 1/6 * (K1laccs + 2*(K2laccs + K3laccs) + K4lacc);
        raccs = raccs + 1/6 * (K1raccs + 2*(K2raccs + K3raccs) + K4racc);
        
        vs_all(i+1,:) = vs;
        ps_all(i+1,:) = ps;
        laccs_all(i+1,:) = laccs;
        raccs_all(i+1,:) = raccs;

        t = t + h;

        if mod(i,10) == 0    
            fprintf("Progress: %i out of %i (%0.3f percent)\n", i, timesteps-1, t/tlast*100)
        end
    end
    
    fprintf("Progress: %i out of %i (100.00 percent)\n", (timesteps-1), (timesteps-1))
end