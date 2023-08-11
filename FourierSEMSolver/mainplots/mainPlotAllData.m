clear all
close all

base_path = '~/data/ARD/simulations';
do_anim = false;

ref = load(sprintf('%s/fourier_L5.00_dx0.010_dt0.000010.mat', '~/data/ARD/simulations'));

sim1 = load(sprintf('%s/dxFourier_div8_dxInterface_div8/fourier_sem_L5.00_dx0.078_dt0.000002_p4_ip6.mat', base_path));
sim2 = load(sprintf('%s/dxFourier_div8_dxInterface_div8/fourier_sem_L5.00_dx0.154_dt0.000004_p4_ip6.mat', base_path));

sims_1 = [sim1,sim2];

sim1 = load(sprintf('%s/dxFourier_x2_dxInterface_x1/fourier_sem_L5.00_dx0.011_dt0.000002_p4_ip6.mat', base_path));
sim2 = load(sprintf('%s/dxFourier_x2_dxInterface_x1/fourier_sem_L5.00_dx0.021_dt0.000006_p4_ip6.mat', base_path));
sim3 = load(sprintf('%s/dxFourier_x2_dxInterface_x1/fourier_sem_L5.00_dx0.043_dt0.000012_p4_ip6.mat', base_path));
sim4 = load(sprintf('%s/dxFourier_x2_dxInterface_x1/fourier_sem_L5.00_dx0.081_dt0.000024_p4_ip6.mat', base_path));
sim5 = load(sprintf('%s/dxFourier_x2_dxInterface_x1/fourier_sem_L5.00_dx0.169_dt0.000049_p4_ip6.mat', base_path));

sims_2 = [sim1,sim2,sim3,sim4,sim5];

sim1 = load(sprintf('%s/dxFourier_x2_dxInterface_x1_dtmax/fourier_sem_L5.00_dx0.011_dt0.000025_p4_ip6.mat', base_path));
sim2 = load(sprintf('%s/dxFourier_x2_dxInterface_x1_dtmax/fourier_sem_L5.00_dx0.021_dt0.000049_p4_ip6.mat', base_path));
sim3 = load(sprintf('%s/dxFourier_x2_dxInterface_x1_dtmax/fourier_sem_L5.00_dx0.043_dt0.000099_p4_ip6.mat', base_path));
sim4 = load(sprintf('%s/dxFourier_x2_dxInterface_x1_dtmax/fourier_sem_L5.00_dx0.081_dt0.000188_p4_ip6.mat', base_path));
sim5 = load(sprintf('%s/dxFourier_x2_dxInterface_x1_dtmax/fourier_sem_L5.00_dx0.169_dt0.000197_p4_ip6.mat', base_path));

sims_3 = [sim1,sim2,sim3,sim4,sim5];

P_order = sim1.scheme_order;
iface_order = sim1.interface_order;

[h1, err1] = calcErr(ref, sims_1);
[h2, err2] = calcErr(ref, sims_2);
[h3, err3] = calcErr(ref, sims_3);

figure(3)
loglog(h2, err2, '-o', 'LineWidth',2)
hold on
loglog(h3, err3, '-v', 'LineWidth',2)
loglog(h1, err1, 'r-o', 'LineWidth',2)
loglog(h2, 1*h2.^2, 'b--', 'LineWidth',2)
loglog(h2(3:end), 10*h2(3:end).^5, 'r--', 'LineWidth',2)
xlabel('\Delta x [m]', 'FontSize', 18)
ylabel('L_2 error', 'FontSize', 18)
title_str = sprintf('FM-SEM_{P=4} convergence');
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

function [h, err] = calcErr(ref, simulations)
    err = zeros(1,length(simulations));
    h = zeros(1,length(simulations));
    calcN = @(p1,p2) length(p1) + length(p2)-1;
    
    nc = @(data) int32(data.tmax/data.dt);
    
    for i=1:length(simulations)        
        sim = simulations(i);
        l1 = sim.x1d_1(end); % same for all simulations
        
        x1d_full1 = [sim.x1d_1(1:end-1), l1+sim.x1d_2];    

        ref1.p = spline(ref.x1d, ref.p, x1d_full1);
        ref1.x1d = x1d_full1; ref1.dt = ref.dt; ref1.tmax = ref.tmax;    

        p1 = sim.p1(nc(sim),:); p2 = sim.p2(nc(sim),:);
        err(i) = norm(abs(ref1.p(nc(ref1),:) - [p1(1:end-1),p2]), 2)/calcN(p1,p2);

        h(i) = sim.dx;
    end
end