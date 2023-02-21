function [z,accs] = rk4_coupled_weq2(dzdt, z0, tmax, timesteps, accs)
    h = tmax/timesteps;
    t = 0.0;
    z = zeros(size(z0,1),timesteps+1);
    z(:,1) = z0;
    
    if nargin == 4
        accs = [];
    end

    for i=1:timesteps
        [K1,K1accs] = dzdt(t,z(:,i),accs);
        [K2,K2accs] = dzdt(t+h/2, z(:,i) + h/2*K1, accs + h/2*K1accs);
        [K3,K3accs] = dzdt(t+h/2, z(:,i) + h/2*K2, accs + h/2*K2accs);
        [K4,K4accs] = dzdt(t+h, z(:,i) + h*K3, accs + h*K3accs);
        z(:,i+1) = z(:,i) + h/6*(K1 + 2*(K2 + K3) + K4);
        accs = accs + h/6*(K1accs + 2*(K2accs + K3accs) + K4accs);
        
        t = t + h;

        if mod(i,10) == 0    
            fprintf("Progress: %i out of %i (%0.3f percent)\n", i, timesteps-1, t/tmax*100)
        end
    end
    
    fprintf("Progress: %i out of %i (100.00 percent)\n", (timesteps-1), (timesteps-1))
end