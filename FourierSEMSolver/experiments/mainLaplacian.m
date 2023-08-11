clear all
close all

P_order = 1;

L = 2*pi;
c = 343;
ppw = 4;
f = 1;
xminmax = [0,L];

dx = 0.1; %c/(ppw*f);

[x1d1,custom1,dxMin1,S, M, Dr] = setupSEMCustom(xminmax, dx, P_order);

%A = [1 3 5; 2 4 6; 7 8 10];

%D2 = Dr'*(A \A)*Dr;
%D2_ = Dr'*Dr;

y = sin(2*pi*f/L*x1d1);
%y = cos(2*pi*f/L*x1d1);
%y = [y(1:end), fliplr(y(1:end-1))];

%y = y(1:2:end);

figure(2)
plot(x1d1, y, 'k-')
hold on
%ddy = 1/(pi)*S*y';
ddy = (custom1.Mx \ custom1.Sx)*y';

factor1 = 1/max(ddy(3:P_order:end-1));
factor2 = 1/max(ddy(2:P_order:end-1));

plot(x1d1(1:end), ddy(1:end), '-or')
%plot(x1d1(1:P_order:end), ddy(1:P_order:end), 'or')
%plot(x1d1(2:P_order:end), -ddy(2:P_order:end), 'oc')
%plot(x1d1(3:P_order:end-1), factor1*ddy(3:P_order:end-1), 'ob', x1d1(2:P_order:end-1), factor2*ddy(2:P_order:end-1), 'ob')
%plot(x1d1(2:P_order:end-1), factor1*ddy(2:P_order:end-1), 'ob')
xlabel('x')
ylabel('sin(2\pi f/L \times x)')
%hl = legend('$\sin$', '$\partial_{xx}\sin$', '$\partial_{xx}\sin$', '$\partial_{xx}\sin$ manually scaled');
hl = legend('$\sin$', '$\partial_{xx}\sin$', '$\partial_{xx}\sin$ manually scaled');
set(hl, 'Interpreter','latex')
set(gca,'FontSize',16)
title(sprintf('$$\\mathrm{D^2} = \\mathrm{M}_x^{-1}\\mathrm{S}_x$$ of order P=%i applied to sine function', P_order), 'FontSize', 12, 'Interpreter','latex')