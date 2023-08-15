clear all
close all

base_path = '~/data/ARD/simulations/dxFourier_div8_dxInterface_div8';
do_anim = false;

ref = load(sprintf('%s/fourier_L5.00_dx0.010_dt0.000010.mat', '~/data/ARD/simulations'));

sim1 = load(sprintf('%s/fourier_sem_L5.00_dx0.077_dt0.000010_p4_ip6.mat', base_path));
sim2 = load(sprintf('%s/fourier_sem_L5.00_dx0.149_dt0.000010_p4_ip6.mat', base_path));

P_order = sim1.scheme_order;
iface_order = sim1.interface_order;

nc = @(data) int32(ref.tmax/data.dt);
l1 = sim1.x1d_1(end); % same for all simulations

x1d_full1 = [sim1.x1d_1(1:end-1), l1+sim1.x1d_2];
x1d_full2 = [sim2.x1d_1(1:end-1), l1+sim2.x1d_2];

ref1.p = spline(ref.x1d, ref.p, x1d_full1);
ref1.x1d = x1d_full1; ref1.dt = ref.dt;
ref2.p = spline(ref.x1d, ref.p, x1d_full2);
ref2.x1d = x1d_full2; ref2.dt = ref.dt;

err = zeros(1,2);

calcN = @(p1,p2) length(p1) + length(p2)-1;

p1 = sim1.p1(nc(sim2),:); p2 = sim1.p2(nc(sim2),:);
err(1) = norm(abs(ref1.p(nc(ref1),:) - [p1(1:end-1),p2]), 2)/calcN(p1,p2);

p1 = sim2.p1(nc(sim2),:); p2 = sim2.p2(nc(sim2),:);
err(2) = norm(abs(ref2.p(nc(ref2),:) - [p1(1:end-1),p2]), 2)/calcN(p1,p2);

h = [sim1.dx,sim2.dx];

figure(3)
loglog(h, err, '-o', 'LineWidth',2)
hold on
loglog(h, 1*h.^5, '--', 'LineWidth',2)
xlabel('\Delta x', 'FontSize', 18)
ylabel('L^2 error', 'FontSize', 18)
title_str = sprintf('FOURIER-SEM convergence \n dx_{fourier} == dx_{sem} / 8  |  dx_{sem-interface} == dx_{sem} / 8');
title(title_str)
%hl = legend(sprintf('SEM P=%i',P_order), '$\mathcal{O}(\Delta x^2)$', '$\mathcal{O}(\Delta x^3)$', '$\mathcal{O}(\Delta x^4)$', 'FontSize', 14);
hl = legend(sprintf('SEM P=%i',P_order), '$\mathcal{O}(\Delta x^5)$', 'FontSize', 14, 'location', 'northwest');
set(hl, 'Interpreter','latex')
set(gca,'FontSize',16)