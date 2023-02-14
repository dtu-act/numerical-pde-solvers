function [psAll] = rk4_weq(dxP, ps0, t0, tlast, timesteps)
    
    psAll = zeros(timesteps, length(ps0));
    psAll(1,:) = ps0;    
    
    ps = zeros(length(ps0),1);
    ps(:) = ps0;    
    
    h = (tlast-t0)/timesteps;
    
    t = 0.0;
    hHalf = h/2.0;
    oneSixth = 1.0/6.0;
    
    for i=1:timesteps-1
        % calculate the pressure and the velocity for the next timestep

        K1p = h * dxP(t,          ps);
        K2p = h * dxP(t + hHalf,  ps + K1p./2.0);      
        K3p = h * dxP(t + hHalf,  ps + K2p./2.0);
        K4p = h * dxP(t + h,      ps + K3p);        
        ps(:) = ps(:) + oneSixth * (K1p + 2*(K2p + K3p) + K4p);
        
        psAll(i+1,:) = ps;
        
        t = t + h;

        if mod(i,10) == 0    
            fprintf("Progress: %i out of %i (%0.3f percent)\n", i, timesteps-1, t/tlast*100)
        end
    end
    
    fprintf("Progress: %i out of %i (100.00 percent)\n", (timesteps-1), (timesteps-1))
end