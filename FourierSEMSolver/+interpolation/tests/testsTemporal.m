%% Test temporal tests resolution adjustment
close all
plotMode = true;

testAdjustTemporal(plotMode);

function testAdjustTemporal(plotMode)
    domain = domainTestSetup(25);    
    x1d = domain.xmin:domain.dx:domain.xmax;
    
    nhistory = 5;     
    
    state = FDTDSimulationTestSetup(x1d, domain, Source(), CustomFDTD(2), nhistory);
    [state_interpl, interpolations] = adjustTemporal(state,domain.dt/2,[1,2],0);

    %testCase = matlab.unittest.TestCase.forInteractiveUse;
    %assertEqual(testCase, p2c_rt, p2c_rt2)
    
    if plotMode
        figure(1)
        tiledlayout(1,1)
        ax1 = nexttile;
        
        hold(ax1, 'on')
        plot(ax1, x1d, state_interpl.p_current, 'o-r');
        plot(ax1, x1d, state.p_prev);
        plot(ax1, x1d, state_interpl.p_prev, 'o-b');
        % nhistory = 5, hence we get 9 interpolations in total
%         plot(ax1, x1d, interpolations(:,1));
%         plot(ax1, x1d, interpolations(:,2));
%         plot(ax1, x1d, interpolations(:,3));
%         plot(ax1, x1d, interpolations(:,4));
%         plot(ax1, x1d, interpolations(:,5));
%         plot(ax1, x1d, interpolations(:,6));
%         plot(ax1, x1d, interpolations(:,7));
%         plot(ax1, x1d, interpolations(:,8));
%         plot(ax1, x1d, interpolations(:,9));
        xlabel(ax1, 'x')
        ylabel(ax1, 'p')
        title(ax1, 'Temporal spline interpolations')
        legend(ax1, 'p_n', 'p_{n-1}', 'Interpolated p_{n-1/2}')
        hold(ax1, 'off')
    end
end