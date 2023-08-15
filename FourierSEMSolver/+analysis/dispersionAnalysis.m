p_fdtd_factor4 = load('../../virtual-acoustics-data/data/p_fdtd_modes_800_l_10_sourcefreq_800_Fn_27440').p_fdtd;
p_fdtd_factor2 = load('../../virtual-acoustics-data/data/p_fdtd_modes_400_l_10_sourcefreq_800_Fn_13720').p_fdtd;

p_wbm_factor2 = load('../../virtual-acoustics-data/data/p_wbm_modes_400_l_10_sourcefreq_800_Fn_13720').p_wbm;
p_wbm_factor4 = load('../../virtual-acoustics-data/data/p_wbm_modes_800_l_10_sourcefreq_800_Fn_27440').p_wbm;

figure(1)
tiledlayout(3,1) % Requires R2019b or later
ax1 = nexttile;
ax2 = nexttile;
ax3 = nexttile;

plot(ax1,p_wbm_factor2); hold(ax1, 'on'); plot(ax1,p_fdtd_factor2); hold off
xlabel(ax1,'n')
ylabel(ax1,'p')
title(ax1,'WBM v.s. FDTD')
ylim(ax1,[-1,1])
legend(ax1,'WBM 1x', 'FDTD 1x')    
hold(ax1, 'off')

plot(ax2,p_fdtd_factor2); hold(ax2, 'on'); plot(ax2,p_fdtd_factor4(1:2:end)); hold off
xlabel(ax2,'n')
ylabel(ax2,'p')
title(ax2,'FDTD comparison for \Delta x resolutions')
ylim(ax2,[-1,1])
legend(ax2,'FDTD 1x', 'FDTD 2x')    
hold(ax2, 'off')

plot(ax3,p_wbm_factor2); hold(ax3, 'on'); plot(ax3,p_wbm_factor4(1:2:end)); hold off
xlabel(ax3,'n')
ylabel(ax3,'p')
title(ax3,'WBM comparison for \Delta x resolutions')
ylim(ax3,[-1,1])
legend(ax3,'WBM 1x', 'WBM 2x')
hold(ax3, 'off')