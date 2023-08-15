function x = idctmodes(X, boundary_type)
    import models.types.*
    
    switch boundary_type
        case BoundaryType.Neumann
            x = idctmodes_neumann(X);
        case BoundaryType.Dirichlet
            x = idctmodes_dirichlet(X);
        otherwise
            error('Boundary condition not supported')
    end
end

function x = idctmodes_neumann(X)
    N = length(X);
    
    Xdouble = [X,fliplr(X(2:end-1))];
    x = real(ifft(Xdouble));
    x = x(1:N)*((N-1)*2);
    
    %fprintf('Im(X) = %f\n', max(imag(X)));
end

function x = idctmodes_dirichlet(X)
    N = length(X);
    
    Xdouble = 1i*[X,-fliplr(X(2:end-1))];
    x = ifft(Xdouble);    
    x = x(1:N)*N;
    
%     fprintf('Im(X) = %f\n', max(imag(Xdouble)));
%     fprintf('Re(X) = %f\n', max(real(Xdouble)));
end

% Using the MATLAB DCT type I, II, III, or IV
function x = idctmodes_matlab(X, type)
    N = length(X);
    
    fprintf('Im(X) = %f\n', max(imag(X)));
    x = idct(X, N, 'Type', type)/(sqrt(2/(N-1)));
end