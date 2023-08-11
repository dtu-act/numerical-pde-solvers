function signal_out = zeroPad(signal, length_out)
    len_zeros = length_out - length(signal);
    signal_out = [signal; zeros(1,len_zeros)'];
end