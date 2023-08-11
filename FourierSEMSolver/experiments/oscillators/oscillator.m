clear all
close all

l = 1;
xs = 0:0.01:l;
n_modes = 20;

for t=0:0.00005:1
    ps = pressure(xs,t,l,n_modes);
    plot(xs,ps)
    ylim([-20 30])
    title("Oscillator spanning 20 modes in 1D domain of length l=1")
    xlabel('length [m]')
    ylabel('Magnitude');
    
    pause(0.1)
end


function [p] = pressure(x,t,l,n_modes)
a=1;
b=1;
c=343;
k = @(i) pi*(i/l);
m = @(i,t) a*exp(1i*c*k(i)*t) + b*exp(-1i*c*k(i)*t); % imaginary part is zero
theta = @(i,x) cos(pi*i/l * x);

p = 0;
for i=0:(n_modes-1)
    p = p + m(i,t)*theta(i,x);
end
end