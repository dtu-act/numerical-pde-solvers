function plotTF(xfreqs,Pfreqs,Pfreqs_ref,param,XY,receiver_pos,sim_type)
    dB = @(amplitude) 20*log10(abs(amplitude));
    
    N = size(receiver_pos,1);
    
    colors = ['b', 'g', 'r', 'm'];
    
    figure()
    if sim_type == "WBM"
        sim_param_str = sprintf('N_w=%i',param);
    else
        sim_param_str = sprintf('P=%i',param);
    end
    
    sg = sgtitle(sprintf('Transfer function, %s, %iHz-%iHz', sim_param_str, min(xfreqs), max(xfreqs)));
    sg.FontSize = 18;
    sg.FontWeight = 'bold';
    
    for k=1:N
        x = receiver_pos(k,1);
        y = receiver_pos(k,2);

        i = meshing.getIndexForMeshCoordinate(XY,x,y);
        
        col = sprintf('%s-', colors(mod(k-1, length(colors))+1));
        
        subplot(N,2,(k-1)*2+1)
        plot(xfreqs,dB(Pfreqs(:,i)), col)
        hold on
        plot(xfreqs,dB(Pfreqs_ref(:,i)), '--r')
        legend(sim_param_str, 'Greens function')
        hold off
        title(sprintf('(r_{x,%i}, r_{y,%i}) = (%0.2f, %0.2f)',i,i,XY(i,1),XY(i,1)))
        xlabel('Frequency [Hz]')
        ylabel('SPL [dB]')
        grid
        set(gca,'fontsize',16)
        
        subplot(N,2,k*2)
        plot(xfreqs,angle(Pfreqs(:,i)), col)
        hold on
        plot(xfreqs,angle(Pfreqs_ref(:,i)), '--r')
        legend(sim_param_str, 'Greens function')
        hold off
        title(sprintf('Angle (Phase), (r_{x,%i}, r_{y,%i}) = (%0.2f, %0.2f)',i,i,XY(i,1),XY(i,1)))
        xlabel('Frequency [Hz]')
        ylabel('Angle (radians)')
        grid
        set(gca,'fontsize',16)
    end
end