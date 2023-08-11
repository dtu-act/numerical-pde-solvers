close all
clear all

[dir_cache, ~, dir_data] = paths.setupPaths("LOCAL");

%==== SIMULATION PARAMETERS =====%
sim_type = "SEM";

f_min_max = [20,1200];
extract_from_to = [20,1200];

lx = 2; ly = 2;
bbox = [0,0; lx,ly];
xy0_rel = [3/5,3/5];

c = 343; rho = 1.225;
aparams = models.AcousticParameters(0,c,rho);

boundary_type = models.BoundaryCondition.Velocity;

if sim_type == "SEM"
    P_order = 4;
    
    NeX = 4;
    NeY = 4;
    xy0 = meshing.calculateSourcePosition(NeX,NeY,bbox,xy0_rel);

    [F_sem_freq, XY, freqs] = caching.loadSimulationData(dir_data, ...
        f_min_max, P_order, NeX, NeY, bbox, xy0_rel, boundary_type, aparams, 'SEM');
    [F_ref_freq, ~, ~]      = caching.loadSimulationData(dir_data, ...
        f_min_max, P_order, NeX, NeY, bbox, xy0_rel, boundary_type, aparams, 'GREENS_REF');
elseif sim_type == "WBM"
    Nwx = 42; Nwy = 42;
    grid_res = 0.1;
    NeX = lx/grid_res; NeY = ly/grid_res;
    
    xy0_rel = [3/5,3/5];
    xy0 = [lx,ly].*xy0_rel;
    
    [F_sem_freq, XY, freqs] = caching.loadSimulationData(dir_data, ...
        f_min_max, 1, NeX, NeY, bbox, xy0_rel, boundary_type, aparams, sprintf('WBM_Nwx%i_Nwy%i', Nwx,Nwy));
    [F_ref_freq, ~, ~]      = caching.loadSimulationData(dir_data, ...
        f_min_max, 1, NeX, NeY, bbox, xy0_rel, boundary_type, aparams, 'GREENS_REF');
else
    error('Simulation type not supported')
end

source_type = "point"; % "point", "gaussian"
mesh_type = "uniform"; % "uniform", "nonuniform"

%=========================%

assert(f_min_max(2) > f_min_max(1))
assert(extract_from_to(2) > extract_from_to(1))

fbin_res = (f_min_max(2)-f_min_max(1))/(length(freqs)-1);
assert(fbin_res == freqs(2)-freqs(1)) % same as above

T = 1/fbin_res;

freq_from = extract_from_to(1);
freq_to = extract_from_to(2);
from_index = find(freqs == freq_from);
to_index = find(freqs == freq_to);
F_sem_freq = F_sem_freq(from_index:to_index,:);
F_ref_freq = F_ref_freq(from_index:to_index,:);
freqs = freqs(from_index:to_index);

%%% PLOT %%
if sim_type == "SEM"
    plotting.plotTF(freqs,F_sem_freq,F_ref_freq,P_order,XY,[0.2,0.2],sim_type)
else
    Nw = Nwx + Nwy;
    plotting.plotTF(freqs,F_sem_freq,F_ref_freq,Nw,XY,[0.2,0.2],sim_type)
end

% pause
% close all
% 
% F_sem_time = validation.frequencyToTimeDomain(F_sem_freq);
% F_ref_time = validation.frequencyToTimeDomain(F_ref_freq);
% plotting.plotTime2D(F_sem_time, F_ref_time, T, XY, [0.2,0.2], 1)