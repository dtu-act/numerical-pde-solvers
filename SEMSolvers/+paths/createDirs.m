function [ref_dir, plot_dir, freq_data_dir] = createDirs(setup_path)    
    
    base_path = setup_path.base_path;
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