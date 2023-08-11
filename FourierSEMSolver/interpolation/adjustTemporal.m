function [state_out,p_t_refined] = adjustTemporal(state_in, dt_fine)
    % return: the current and previous time step interpolated at dt_fine
    
    dt_coarse = state_in.domain.dt;
    
    assert(rem(dt_coarse, dt_fine) == 0)
    
    if dt_coarse <= dt_fine
        error("error: dt_coarse < dt_fine")
    end

    p_history = state_in.p_history;

    steps_in = size(p_history,2);
    steps_out = 2;

    x = 0:dt_coarse:(steps_in-1)*dt_coarse; % number of steps of timesteps (history)
    xx = 0:dt_fine:(steps_out-1)*dt_fine;
    
    p_t_refined = spline(x,p_history,xx);

    state_out = copy(state_in);
    state_out.domain.dt = dt_fine;
    state_out.adjust(p_t_refined(:,1)', p_t_refined(:,2)');
end