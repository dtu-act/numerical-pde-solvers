function [z] = rk4_coupled_weq2(dzdt, z0, tmax, timesteps)   
    h = tmax/timesteps;
    t = 0.0;
    z = zeros(size(z0,1),timesteps+1);
    z(:,1) = z0;

    for i=1:timesteps
        K1 = dzdt(t,z(:,i));
        K2 = dzdt(t+h/2, z(:,i) + h/2*K1);
        K3 = dzdt(t+h/2, z(:,i) + h/2*K2);
        K4 = dzdt(t+h, z(:,i) + h*K3);
        z(:,i+1) = z(:,i) + h/6*(K1 + 2*(K2 + K3) + K4);
        
        t = t + h;

        if mod(i,10) == 0    
            fprintf("Progress: %i out of %i (%0.3f percent)\n", i, timesteps-1, t/tmax*100)
        end
    end
    
    fprintf("Progress: %i out of %i (100.00 percent)\n", (timesteps-1), (timesteps-1))
end