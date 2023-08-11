function [ref_dir, plot_dir, freq_data_dir] = setupPaths(env)
    if env == "LOCAL"
        setup_paths = jsondecode(fileread('../scripts/simulation_setups/simulation_setup_local.json'));
    elseif env == "HPC"
        addpath('~/matlab/distmesh-master') % use this running on HPC
        setup_paths = jsondecode(fileread('../scripts/simulation_setups/simulation_setup_hpc.json'));
    else
        error('Environment not supported')
    end
    
    base_path = setup_paths.base_path;
    plot_dir = sprintf('%s/plots', base_path);
    ref_dir = sprintf('%s/ref',base_path);
    freq_data_dir = sprintf('%s/freq_data',base_path);
    
    if ~exist(plot_dir, 'dir')
        mkdir(plot_dir)
    end
    
    if ~exist(ref_dir, 'dir')
        mkdir(ref_dir)
    end
    
    if ~exist(freq_data_dir, 'dir')
        mkdir(freq_data_dir)
    end
end