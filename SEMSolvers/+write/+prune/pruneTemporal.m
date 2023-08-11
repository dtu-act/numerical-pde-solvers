function [tsteps_out,p_out,dt_out] = pruneTemporal(dt_in,fmax,tsteps_in,p_in,ppw_out)
    dt_out = 1.0/(fmax*ppw_out);
    tsteps_write = max(floor(dt_out/dt_in), 1);

    tsteps_out = tsteps_in(1:tsteps_write:end);
    p_out = p_in(:,:,1:tsteps_write:end);
    dt_out = tsteps_out(2)-tsteps_out(1);
end