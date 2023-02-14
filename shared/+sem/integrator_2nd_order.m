function [ps_all] = integrator_2nd_order(dxP, ps0_current, ps0_prev, dt, timesteps, Q)
    
    ps_all = zeros(timesteps, length(ps0_current));
    
    p_current = zeros(length(ps0_current),1);
    p_prev = zeros(length(ps0_prev),1);
    
    p_current(:) = ps0_current;
    p_prev(:) = ps0_prev;
    
    ps_all(1,:) = ps0_prev;
    
    for i=1:timesteps-1
        % calculate the pressure and the velocity for the next timestep        
        ps_all(i+1,:) = 2*p_current - p_prev + dt^2 * dxP(i*dt, p_current) + Q(i*dt);
        p_prev = p_current;
        p_current(:) = ps_all(i+1,:);

        if mod(i,10) == 0    
            fprintf("Progress: %i out of %i (%0.3f percent)\n", i, timesteps-1, i/timesteps*100)
        end
    end
    
    fprintf("Progress: %i out of %i (100.00 percent)\n", (timesteps-1), (timesteps-1))
end