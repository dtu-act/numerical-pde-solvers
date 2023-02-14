function [Ltilde,btilde] = assemblyImpedanceBounds2D(L,b,x,y,conn,aparams,bdetect_f,impedance_f)
    Ltilde = L; btilde = b;

    bound_indexes = [];
    
    indxs = unique(conn);
    for i=indxs'
        if bdetect_f(x(i), y(i)) == 0
            bound_indexes = [bound_indexes i]; % TODO: slow implementation
        end        
    end

    for i=bound_indexes
        % normalized impedance? 1i*omega*rho = 1 
        btilde(:) = btilde(:) - 1i*aparams.omega*aparams.rho/impedance_f(x(i), y(i)) * Ltilde(:,i);
        Ltilde(i,:) = 0;
        Ltilde(:,i) = 0;
        Ltilde(i,i) = 1;
    end
    
    for i=bound_indexes
        %fprintf('(x,y)=(%0.3f,%0.3f) f(x,y)=%0.3f\n',x(i),y(i), bound_cond_f(x(i),y(i)));
        btilde(i) = impedance_f(x(i),y(i));
    end
end