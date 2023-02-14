function [Ltilde,btilde] = assemblyDirichletBounds2D(L,b,x,y,conn,bdetect_f,pressure_f)
    Ltilde = L; btilde = b;

    bound_indexes = [];
    
    indxs = unique(conn);
    for i=indxs'
        if bdetect_f(x(i), y(i)) == 0
            bound_indexes = [bound_indexes i]; % TODO: slow implementation
        end        
    end

    for i=bound_indexes
        btilde(:) = btilde(:) - pressure_f(x(i), y(i)) * Ltilde(:,i);
        Ltilde(i,:) = 0;
        Ltilde(:,i) = 0;
        Ltilde(i,i) = 1;
    end
    
    for i=bound_indexes
        %fprintf('(x,y)=(%0.3f,%0.3f) f(x,y)=%0.3f\n',x(i),y(i), bound_cond_f(x(i),y(i)));
        btilde(i) = pressure_f(x(i),y(i));
    end
end