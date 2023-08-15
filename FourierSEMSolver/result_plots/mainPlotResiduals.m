clear all
close all

s1 = load('~/data/fm-sem/residual_removed_order1_1_order2_2_n500.mat');
s2 = load('~/data/fm-sem/residual_removed_order1_2_order2_2_n500.mat');

x1d_ref = s1.x1d_2;
x1d_1 = s1.x1d_1;
x1d_2 = s2.x1d_1;

p_ref = s1.p2;
p1 = s1.p1;
p2 = s2.p1;

dxxP_ref = s1.dxxP2;
dxxP1 = s1.dxxP1;
dxxP2 = s2.dxxP1;

h1 = figure(1);
set(h1,'Position',[10 10 500 300])
h2 = figure(2);
set(h2,'Position',[10 10 500 300])

l1 = 'FDTD reference, 2nd order';
l2 = 'SEM P=1';
l3 = 'SEM P=2';

figure(1)
plot(x1d_ref,p_ref,'--','LineWidth', 2)
hold on
plot(x1d_1,p1,'-or','LineWidth', 2)
plot(x1d_2,p2,'-ob','LineWidth', 2)
hold off

%title('Left boundary')
legend(l1,l2,l3);
xlabel('x [m]')
ylabel('Reflected pressure [Pa]')
ylim([0,0.6])
%xlim([min(x1d_1(1),x1d_2(1)),max(x1d_1(end),x1d_2(end))])
xlim([min(x1d_1(1),x1d_2(1)),1])
title('Left boundary vicinity')
set(gca,'FontSize',16)

figure(2)
plot(x1d_ref(2:end),dxxP_ref,'--','LineWidth', 2)
hold on
plot(x1d_1(2:end),dxxP1,'-or','LineWidth', 2)
plot(x1d_2(2:end),dxxP2,'-ob','LineWidth', 2)
hold off                
legend(l1,l2,l3,'location', 'southeast');
xlabel('x [m]')
title('Laplacian with Neumann boundaries')
%title('Laplacians', 'fontsize', 16)
xlim([0,0.4])
set(gca,'FontSize',16)