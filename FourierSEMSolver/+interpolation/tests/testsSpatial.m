clear all
close all

%% Test temporal resolution adjustment
figure(1)
tiledlayout(2,1)
ax1 = nexttile;
ax2 = nexttile;
    
testFinerSpatialInterpolation(ax1);
testCoarserSpatialInterpolation(ax2);

function testFinerSpatialInterpolation(ax)
    custom = CustomFDTD(2);
    
    domain_orig = domainTestSetup(50);     
    x1d = domain_orig.xmin:domain_orig.dx:domain_orig.xmax;
    state = FDTDSimulationTestSetup(x1d, domain_orig, zeroSource(x1d), custom, 5);
    
    dx_i = domain_orig.dx*1/2;
    state_i = adjustSpatial(state, dx_i);    
    x1d_i = domain_orig.xmin:dx_i:domain_orig.xmax;
    
    % compare
    domain_compare = domainTestSetup(100);
    state_compare = FDTDSimulationTestSetup(x1d_i, domain_compare, zeroSource(x1d_i), custom, 5);
    
    plotStates(x1d, x1d_i, state, state_i, state_compare, ax, 'Refined spatial spline interpolations')
end

function testCoarserSpatialInterpolation(ax)
    custom = CustomFDTD(2);
    
    domain_orig = domainTestSetup(50);
    x1d = domain_orig.xmin:domain_orig.dx:domain_orig.xmax;
    sim = FDTDSimulationTestSetup(x1d, domain_orig, zeroSource(x1d), custom, 5);
    
    dx_i = domain_orig.dx*2;
    state_i = adjustSpatial(sim, dx_i);
    x1d_i = domain_orig.xmin:dx_i:domain_orig.xmax;
    
    % compare
    domain_exact = domainTestSetup(25);
    sim_exact = FDTDSimulationTestSetup(x1d_i, domain_exact, zeroSource(x1d_i), custom, 5);
    
    plotStates(x1d, x1d_i, sim, state_i, sim_exact, ax, 'Coarsened spatial spline interpolations')
end

function plotStates(x1d, x1d_i, state, state_i, state_compare, ax, title_text)
    plot(ax, x1d, state.p_current, 'or');
    hold(ax, 'on')
    plot(ax, x1d_i, state_i.p_current, '--xb');
    plot(ax, x1d_i, state_compare.p_current, '-r');
    plot(ax, x1d, state.p_prev, 'og');
    plot(ax, x1d_i, state_i.p_prev, '--xk');
    plot(ax, x1d_i, state_compare.p_prev, '-g');
    xlabel(ax, 'x')
    ylabel(ax, 'p')
    title(ax, title_text)
    legend(ax, 'p^{n+1} \Deltax', 'p^{n+1} \Deltax_{interp}', 'p^{n+1} \Deltax_{exact}', .../
        'p^{n} \Deltax', 'p^{n} \Deltax_{interp}', 'p^{n} \Deltax_{exact}')
    hold(ax, 'off')        
end