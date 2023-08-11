function [wavenumber, scf, wavecount] = wavenumber_gen2D(domain, params, wx_max, wy_max)

    k = params.k;

    lx = domain.lx;
    ly = domain.ly;

    % WAVE NUMBER GENERATION
    Kxwr = [];
    Kywr = [];
    Kxws = [];
    Kyws = [];

    for w1 = 0:wy_max
        Ky = w1*pi/ly;
        Kx = sqrt(k^2-Ky^2); % can be complex
        Kxwr = [ Kxwr Kx -Kx ];
        Kywr = [ Kywr Ky  Ky ];
    end
    
    for w2 = 0:wx_max
        Kx=w2*pi/lx;
        Ky=sqrt(k^2-Kx^2); % can be complex
        Kxws = [ Kxws Kx  Kx ];
        Kyws = [ Kyws Ky -Ky ];
    end
    
    % number of wave function at each subdomain
    nwr = length(Kxwr);
    nws = length(Kxws);
    nw = nwr + nws;

    %  for the scale factor in wave function 
    scf_r = zeros(1,nwr);
    scf_s = zeros(1,nws);

    for r = 1:nwr
        if imag(Kxwr(r)) > 0
            scf_r(r) = 1;
        end
    end

    for s = 1:nws
        if imag(Kyws(s)) > 0
            scf_s(s) = 1;
        end
    end

    wavenumber = struct('Kxwr',Kxwr, 'Kywr',Kywr, 'Kxws',Kxws, 'Kyws',Kyws);
    scf        = struct('scf_r',scf_r, 'scf_s',scf_s);
    wavecount  = struct('nwr',nwr, 'nws',nws, 'nw',nw);
end