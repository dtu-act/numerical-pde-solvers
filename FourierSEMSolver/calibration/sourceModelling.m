clear all

c = 343;

Fs = 1750;
fmax_pulse = 100;
plot_from = 20;
dt = 1/Fs;
dct_type = 2;
boundary = BoundaryType.Neumann;

[source_t, source_x, mu] = GaussianSourceDiff(dt,0,c,fmax_pulse);

ns = 0:1000;
s = source_t(ns);

plotSignalPropertiesDD(ns*dt, source_t(ns), boundary, Fs, fmax_pulse*2, plot_from, dct_type, 1);