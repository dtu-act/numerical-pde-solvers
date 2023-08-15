function X = dctmodes(x, boundary_type)
    import models.types.*
    
    switch boundary_type
        case BoundaryType.Neumann
            X = dctmodes_neumann(x);
        case BoundaryType.Dirichlet
            X = dctmodes_dirichlet(x);
        otherwise
            error('Boundary condition not supported')
    end
end

function X = dctmodes_neumann(x)    
    xdouble = [x, fliplr(x(2:end-1))];
    N = length(xdouble);

    Xdouble = fft(xdouble);
    X = real(Xdouble(1:length(x)));
    X = X/N;
    
    %fprintf('Im(X) = %f\n', max(imag(Xdouble)));
end

function X = dctmodes_dirichlet(x)    
    N = length(x);
    xdouble = [x, -fliplr(x(2:end-1))];
    
    Xdouble = fft(xdouble);
    X = imag(Xdouble(1:N));
    X = X/N;
    
%     fprintf('Im(X) = %f\n', max(imag(Xdouble)));
%     fprintf('Re(X) = %f\n', max(real(Xdouble)));
end

% Using the MATLAB DCT type I, II, III, or IV
function X = dctmodes_matlab(x, type)
    N = length(x);
    X = dct(x, N, 'Type', type)*(sqrt(2/(N-1)));
end

% https://dsp.stackexchange.com/questions/49632/dft-fft-of-a-real-even-function-doesnt-yield-real-only-dft-signal