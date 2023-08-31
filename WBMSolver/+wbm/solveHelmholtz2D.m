function [P, xs, ys, wavecount] = solveHelmholtz2D(domains, acoustic_params, sim_params, grid_res)
    
    if nargin == 3
        grid_res = 0.1;
    end
    
    N_domains = length(domains);
    
    return_as_cell = true;
    
    if N_domains == 1 && ~iscell(domains)
        return_as_cell = false;
        domains = {domains};
    end
    
    wavenumbers = cell(N_domains,1);
    scfs = cell(N_domains,1);
    wavecounts = cell(N_domains,1);
    
    trunc_par = sim_params.truncation_param;
    
    % Calculate wavefunctions
    for i=1:N_domains
        domain = domains{i};
        
        % Truncation of the wave function
        wx_max = max(2, ceil( domain.lx * trunc_par * acoustic_params.k/pi )); % number of wave function r of each domain
        wy_max = max(2, ceil( domain.ly * trunc_par * acoustic_params.k/pi )); % number of wave function s of each domain

        [wavenumber, scf, wavecount] = wbm.wavenumber_gen2D(domain, acoustic_params, wx_max, wy_max);        
        wavenumbers{i} = wavenumber; scfs{i} = scf; wavecounts{i} = wavecount;
    end
    
    % Initialize the system matrices
    A = cell(N_domains, N_domains);
    b = cell(N_domains, 1);
    
    for i=1:N_domains
        for j=1:N_domains
            A{i,j} = zeros(wavecounts{i}.nw, wavecounts{j}.nw);
        end
        b{i} = zeros(wavecounts{i}.nw,1);
    end

    for i=1:N_domains
        domain = domains{i};
        fprintf('D_%i\n',i)
        
        wavenumber = wavenumbers{i}; scf = scfs{i}; wavecount = wavecounts{i};
        
        disp('Boundaries')
        [Ai, Fi] = wbm.assembly.assembly2D(domain, wavenumber, scf, wavecount, trunc_par, acoustic_params);
        
        disp('Interface')
        [Ai_I, Fi_I] = wbm.assembly.assemblyInterface2D(domain, wavenumber, scf, wavecount, trunc_par, acoustic_params);
        
        A{i,i} = Ai + Ai_I;
        b{i} = Fi + Fi_I;
    end
    
    % we assume that all domains share an interface
    for i=1:N_domains        
        for j=1:N_domains
            if i ~= j
                fprintf('D_{%i,%i}\n',i,j)
                disp('Coupling')

                [Cij, Fij] = wbm.assembly.assemblyCoupling2D(i,j,domains,wavenumbers,scfs,wavecounts,trunc_par,acoustic_params);

                A{i,j} = Cij;
                b{i} = b{i} + Fij;
            end
        end
    end
    
    AA = cell2mat(A);
    bb = cell2mat(b);

    Pw = AA \ bb;
    
    i_to = 0;
    
    P = cell(1,N_domains); xs = cell(1,N_domains); ys = cell(1,N_domains);
    for i=1:N_domains
        wavenumber = wavenumbers{i}; scf = scfs{i}; wavecount = wavecounts{i};
        
        i_from = i_to + 1;
        i_to = i_from + wavecount.nw - 1;
        
        [Pi, xsi, ysi] = wbm.postprocess2D(Pw(i_from:i_to), domains{i}, acoustic_params, wavenumber, scf, grid_res);
        
        P{i} = Pi;
        xs{i} = xsi;
        ys{i} = ysi;
    end
    
    assert(i_to == length(Pw))
    
    if ~return_as_cell
        P = P{1};
        xs = xs{1};
        ys = ys{1};
    end
end