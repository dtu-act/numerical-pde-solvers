close all
clear all

% PLOT: CPU w.r.t. domain size
x = [6,13,25,50];
y100_ref = [0.7,0.9,1.4,2.1];
y100_sim = [2.1,6.6,29.6,149.2];

% PLOT: CPU w.r.t. domain size
x_frac = [10,20,50,70,80,90,95];
yfrac_ref = [149.2,149.2,149.2,149.2,149.2,149.2,149.2];
yfrac_sim = [185.9,148.390,53.310,21.147,15.426,10.745,8.385];

figure(1)
loglog(x, y100_ref, '-ob', 'LineWidth',2)
hold on
loglog(x, y100_sim, '-or', 'LineWidth',2)
loglog(x, 0.1*x.*log(x), 'b--', 'LineWidth',2)
loglog(x, 0.1*x.^2, 'r--', 'LineWidth',2)
for i=1:length(x)
    text(x(i),y100_ref(i)*1.3,sprintf('%ix',round(y100_sim(i)/y100_ref(i))), 'FontSize', 14);
end
xlabel('Length (l) [m]', 'FontSize', 18)
ylabel('Time [s]', 'FontSize', 18)
title_str = sprintf('CPU time: FM vs. SEM in full domain');
title(title_str)
hl = legend('FM 100\%', 'SEM 100\%', '$\mathcal{O}(L\log(L))$', '$\mathcal{O}(L^2)$');
set(hl, 'Interpreter','latex')
set(gca,'FontSize',16)
hold off

figure(2)
loglog(x_frac, yfrac_sim, '-ro', 'LineWidth',2)
hold on
loglog(x_frac, yfrac_ref, '--b', 'LineWidth',2)
loglog([50,60,80,95], 5*10e5*1./[50,60,80,95].^3, 'r--', 'LineWidth',2)
for i=1:length(x_frac)
    if i==1
        offset_x = 1.0;
        offset_y = 0.9;
    elseif i==2
        offset_x = 0.95;
        offset_y = 0.85;
    else
        offset_x = 1.0;
        offset_y = 1.2;
    end
    text(offset_x*x_frac(i),yfrac_sim(i)*offset_y,sprintf('%0.1fx',yfrac_ref(i)/yfrac_sim(i)), 'FontSize', 14);
end
xlabel('FM partition size (r) [%]', 'FontSize', 18)
ylabel('Time [s]', 'FontSize', 18)
title_str = sprintf('CPU time: FM-SEM coupling in L=50 m domain');
title(title_str)
hl = legend('FM-SEM coupling', 'SEM full domain', '$\mathcal{O}(r^3)$', 'location', 'southwest');
set(hl, 'Interpreter','latex')
set(gca,'FontSize',16)
hold off
return

loglog(h2, err2, '-o', 'LineWidth',2)
hold on
loglog(h3, err3, '-v', 'LineWidth',2)
loglog(h1, err1, 'r-o', 'LineWidth',2)
loglog(h2, 1*h2.^2, 'b--', 'LineWidth',2)
loglog(h2(3:end), 10*h2(3:end).^5, 'r--', 'LineWidth',2)
xlabel('\Delta x', 'FontSize', 18)
ylabel('L^2 error', 'FontSize', 18)
title_str = sprintf('FOURIER-SEM P=4 convergence');
title(title_str)
% hl = legend('$1/2\Delta x_{Fourier} = \Delta x_{SEM}$', ...\
%     '$1/2\Delta x_{Fourier} = \Delta x_{SEM}, 1/2\Delta t_{Fourier} = \Delta t_{SEM}$', ...\
%     '$\Delta x_{Fourier} = \Delta x_{SEM-interface} = 8\Delta x_{SEM}$', ...\
%     '$\mathcal{O}(\Delta x^2)$', '$\mathcal{O}(\Delta x^5)$', 'FontSize', 14, 'location', 'northwest');
hl = legend('$\Delta x$ 1:2, $\Delta t$ 1:1', ...\
    '$\Delta x$ 1:2, $\Delta t$ 2:1', ...\
    'SEM P=1 layer $\Delta x$ refined $x8$', ...\
    '$\mathcal{O}(\Delta x^2)$', '$\mathcal{O}(\Delta x^5)$', 'FontSize', 14, 'location', 'northwest');
set(hl, 'Interpreter','latex')
set(gca,'FontSize',16)