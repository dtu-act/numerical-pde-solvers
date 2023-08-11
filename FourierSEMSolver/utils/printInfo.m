function printInfo(state, s_t, method, fig)
    params = state.params;
    params_s = state.params_source;
    
    fprintf("-------------------\n");
    fprintf("%s\n", method);
    fprintf("-------------------\n");
    fprintf("Length = %f [m]\n", params.l);
    fprintf("k = %i [modes]\n", params.nmodes);
    fprintf("dt = %f \n", params.dt);
    fprintf("dx = %f [m]\n", params.dx);
    if ~isempty(params_s)
        fprintf("Nyquist frequency source: %f [Hz]\n", params_s.fn_source);          
    end
    fprintf("Nyquist frequency: %f [Hz]\n\n", params.fn);
    
    if ~isempty(params_s)
        plot_source_from = 20;  % [Hz]
        plotSignalPropertiesDD((1:100)*params.dt, s_t(1:100), state.domain.boundary, params.fs, params_s.fn_source*2, ...
                               plot_source_from, 1, fig)
    end    
end