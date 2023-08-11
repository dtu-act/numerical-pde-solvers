clear all
close all

base_path = '~/data/ARD/simulations/dxFourier_x2_dxInterface_x1_dtmax';
do_anim = false;

ref = load(sprintf('%s/fourier_L5.00_dx0.010_dt0.000010.mat', '~/data/ARD/simulations'));

sim1 = load(sprintf('%s/fourier_sem_L5.00_dx0.005_dt0.000012_p4_ip6.mat', base_path));
sim2 = load(sprintf('%s/fourier_sem_L5.00_dx0.011_dt0.000025_p4_ip6.mat', base_path));
sim3 = load(sprintf('%s/fourier_sem_L5.00_dx0.021_dt0.000049_p4_ip6.mat', base_path));
sim4 = load(sprintf('%s/fourier_sem_L5.00_dx0.043_dt0.000098_p4_ip6.mat', base_path));
sim5 = load(sprintf('%s/fourier_sem_L5.00_dx0.081_dt0.000185_p4_ip6.mat', base_path));

P_order = sim1.scheme_order;
iface_order = sim1.interface_order;

nc = @(data) int32(data.tmax/data.dt);

l1 = sim1.x1d_1(end); % same for all simulations

x1d_full1 = [sim1.x1d_1(1:end-1), l1+sim1.x1d_2];
x1d_full2 = [sim2.x1d_1(1:end-1), l1+sim2.x1d_2];
x1d_full3 = [sim3.x1d_1(1:end-1), l1+sim3.x1d_2];
x1d_full4 = [sim4.x1d_1(1:end-1), l1+sim4.x1d_2];
x1d_full5 = [sim5.x1d_1(1:end-1), l1+sim5.x1d_2];

ref1.p = spline(ref.x1d, ref.p, x1d_full1);
ref1.x1d = x1d_full1; ref1.dt = ref.dt; ref1.tmax = ref.tmax;
ref2.p = spline(ref.x1d, ref.p, x1d_full2);
ref2.x1d = x1d_full2; ref2.dt = ref.dt; ref2.tmax = ref.tmax;
ref3.p = spline(ref.x1d, ref.p, x1d_full3);
ref3.x1d = x1d_full3; ref3.dt = ref.dt; ref3.tmax = ref.tmax;
ref4.p = spline(ref.x1d, ref.p, x1d_full4);
ref4.x1d = x1d_full4; ref4.dt = ref.dt; ref4.tmax = ref.tmax;
ref5.p = spline(ref.x1d, ref.p, x1d_full5);
ref5.x1d = x1d_full5; ref5.dt = ref.dt; ref5.tmax = ref.tmax;

sim = sim1;
if do_anim
    figure(5)
    ref_t = spline(ref.x1d, ref.p, x1d_full5);
    t_factor = ref.dt/sim.dt;
    for n=1:nc(sim)
        if mod(n,5) == 0 || n == iter
            plot(ref.x1d, ref.p(n,:), '--')
            hold on
            plot([sim.x1d_1(1:end-1), l1+sim.x1d_2], [sim.p1(n,1:end-1), sim.p2(n,:)], 'o-')             
            hold off

            xlabel('x')
            ylabel('p')
            ylim([-1,1])
            drawnow
        end
    end
end

%% PLOT
figure(1)
plot(ref.x1d, ref.p(nc(ref),:), 'r--')
hold on
plot([sim1.x1d_1(1:end-1), l1+sim1.x1d_2], [sim1.p1(nc(sim1),1:end-1), sim1.p2(nc(sim1),:)], '-')
plot([sim2.x1d_1(1:end-1), l1+sim2.x1d_2], [sim2.p1(nc(sim2),1:end-1), sim2.p2(nc(sim2),:)], '-')
plot([sim3.x1d_1(1:end-1), l1+sim3.x1d_2], [sim3.p1(nc(sim3),1:end-1), sim3.p2(nc(sim3),:)], '-')
plot([sim4.x1d_1(1:end-1), l1+sim4.x1d_2], [sim4.p1(nc(sim4),1:end-1), sim4.p2(nc(sim4),:)], '-')
plot([sim5.x1d_1(1:end-1), l1+sim5.x1d_2], [sim5.p1(nc(sim5),1:end-1), sim5.p2(nc(sim5),:)], '-')
xline(l1, '--', 'LineWidth',2)
hold off
legend('Reference', ...    
    sprintf('SEM dx=%0.3f',sim1.dx), ...\
    sprintf('SEM dx=%0.3f',sim2.dx), ...\
    sprintf('SEM dx=%0.3f',sim3.dx), ...\
    sprintf('SEM dx=%0.3f',sim4.dx), ...\
    sprintf('SEM dx=%0.3f',sim5.dx))
title(sprintf('Convergence at time t=%0.3f', ref.tmax))
xlabel('x')
ylabel('pressure')

%% PLOT
sim = sim1;
figure(2)
plot(ref.x1d, ref.p(nc(ref),:), 'r--')
hold on
plot([sim.x1d_1(1:end-1), l1+sim.x1d_2], [sim.p1(nc(sim),1:end-1), sim.p2(nc(sim),:)], 'bo')
xline(l1, '--', 'LineWidth',2)
hold off
legend('Reference', sprintf('SEM dx=%0.3f',sim.dx))
title(sprintf('Pressure at time t=%0.3f', ref.tmax))
xlabel('x')
ylabel('pressure')

err = zeros(1,3);

calcN = @(p1,p2) length(p1) + length(p2)-1;

p1 = sim1.p1(nc(sim1),:); p2 = sim1.p2(nc(sim1),:);
err(1) = norm(abs(ref1.p(nc(ref1),:) - [p1(1:end-1),p2]), 2)/calcN(p1,p2);

p1 = sim2.p1(nc(sim2),:); p2 = sim2.p2(nc(sim2),:);
err(2) = norm(abs(ref2.p(nc(ref2),:) - [p1(1:end-1),p2]), 2)/calcN(p1,p2);

p1 = sim3.p1(nc(sim3),:); p2 = sim3.p2(nc(sim3),:);
err(3) = norm(abs(ref3.p(nc(ref3),:) - [p1(1:end-1),p2]), 2)/calcN(p1,p2);

p1 = sim4.p1(nc(sim4),:); p2 = sim4.p2(nc(sim4),:);
err(4) = norm(abs(ref4.p(nc(ref4),:) - [p1(1:end-1),p2]), 2)/calcN(p1,p2);

p1 = sim5.p1(nc(sim5),:); p2 = sim5.p2(nc(sim5),:);
err(5) = norm(abs(ref5.p(nc(ref5),:) - [p1(1:end-1),p2]), 2)/calcN(p1,p2);

h = [sim1.dx,sim2.dx,sim3.dx,sim4.dx,sim5.dx];

save(sprintf('%s/convergence_data_x1_x1.mat', base_path), 'h', 'err')

figure(3)
loglog(h, err, '-o', 'LineWidth',2)
hold on
%loglog(h, 1*h.^2, '--', 'LineWidth',2)
%loglog(h, 0.1*h.^3, '--', 'LineWidth',2)
loglog(h, 10000*h.^5, '--', 'LineWidth',2)
xlabel('\Delta x', 'FontSize', 18)
ylabel('L^2 error', 'FontSize', 18)
title_str = sprintf('FOURIER-SEM convergence \n dx_{fourier} == dx_{sem} x 2  |  dx_{sem-interface} == dx_{sem}/2');
title(title_str)
%hl = legend(sprintf('SEM P=%i',P_order), '$\mathcal{O}(\Delta x^2)$', '$\mathcal{O}(\Delta x^3)$', '$\mathcal{O}(\Delta x^4)$', 'FontSize', 14);
hl = legend(sprintf('SEM P=%i',P_order), '$\mathcal{O}(\Delta x^2)$', 'FontSize', 14, 'location', 'northwest');
set(hl, 'Interpreter','latex')
set(gca,'FontSize',16)