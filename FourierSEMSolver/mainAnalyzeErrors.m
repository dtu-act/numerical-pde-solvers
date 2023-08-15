clear all
close all

import utilsDD.*
import plotting.*

base_path = '~/data/fm-sem/';
plots_path = sprintf('%s/%s', base_path, 'plots');
mkdir(base_path)
mkdir(plots_path)

write_gif = false;
write_plot = true && ~write_gif;
if write_gif
    overwrite_time_at = -1;
else
    overwrite_time_at = 0.19;
end

%data = load(sprintf('%s/FOURIER_SEM_ppw4_8_src_LEFT_fmax_1000.mat', base_path));
data = load(sprintf('%s/FOURIER_SEM_ppw4_4_src_LEFT_fmax_1000.mat', base_path));

iter = data.iter;
fmax = data.fmax;
c = 343;

p1 = data.p1; p2 = data.p2; p_ref = data.p_ref;
s1 = data.s1; s2 = data.s2; s_ref = data.s_ref;
x1d_1 = s1.domain.x1d; x1d_ref = s_ref.domain.x1d;
x1d_2 = s1.domain.x1d(end) - s2.domain.x1d(1) + s2.domain.x1d;

p_sim = spline([x1d_1(1:end-1),x1d_2], [p1(:,1:end-1),p2(:,:)], x1d_ref);

tnorm = 0.007; % elapsed time to measure pressure used for normalization
n_norm = round(tnorm/s1.domain.dt);
p_sim = 0.5*p_sim/max(p_sim(n_norm,:));
p_ref = 0.5*p_ref/max(p_ref(n_norm,:));

tmax = iter*s1.domain.dt;
assert(abs(tmax - iter*s_ref.domain.dt) < eps('double'))

% find receiver position
r0_pos = 6.0;
rec_i = getClosestIndex(x1d_ref, r0_pos);

% RMS relative error
sim_rms = rms(p_sim(:,rec_i));
ref_rms = rms(p_ref(:,rec_i));
rel_rms_err = abs(sim_rms - ref_rms)/ref_rms;

% relative error
db_threshold = -60;
[i1_db] = find(20*log(abs(p_ref(:,rec_i))) > db_threshold);
[i2_db] = find(20*log(abs(p_sim(:,rec_i))) > db_threshold);
        
rel_err1 = mean(abs(p_ref(i1_db,rec_i) - p_sim(i1_db,rec_i))./abs(p_ref(i1_db,rec_i)));
%rel_err2 = mean(abs(p_ref(i2_db,rec_i) - p_sim(i2_db,rec_i))./abs(p_ref(i2_db,rec_i)));
rel_err2 = rel_err1;

%rel_err2 = abs(p_ref(i2_db,rec_i) - p_sim(i2_db,rec_i))./abs(p_ref(i2_db,rec_i));
%rel_err2(397) = 0;
%rel_err2(472) = 0;
%rel_err2 = mean(rel_err2);

% max error
max_err = max(abs(p_ref(:,rec_i) - p_sim(:,rec_i)));

fprintf('rms relative err: %0.4f\n', rel_rms_err*100)
fprintf('relative err: %0.4f \n', max(rel_err1,rel_err2)*100)
fprintf('max err: %0.4f\n',max_err)

%% PLOT
h = figure(1);
set(h,'Position',[10 10 1000 500])
if write_gif
    % only 3 tiles, we do not plot frequency plot
    tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'none')
    ax1 = nexttile;
    ax2 = nexttile;
    ax3 = nexttile;
else
    tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'none')
    ax1 = nexttile;
    ax2 = nexttile;
    ax3 = nexttile;
    ax4 = nexttile;
end


tvec = (1:iter)*s1.domain.dt;
[source, src_x, src_t] = sources.gaussianSourceDiff1D(c, fmax);

if write_gif
    from_iter = 1;
    max_iter = iter;
else
    from_iter = iter;
    max_iter = from_iter;
end

for n=from_iter:iter
    t = n*s1.domain.dt;
    
    if  overwrite_time_at > 0
        n_i = round(overwrite_time_at/s1.domain.dt);
    else
        n_i = n;
    end
    
    if mod(n-1,10) == 0 || n==max_iter
%         i1 = find(20*log(abs(p_ref(n,:))) > db_threshold);
%         i2 = find(20*log(abs(p_sim(n,:))) > db_threshold);
%         plot(ax1,x1d_ref(i1),p_sim(n,i1),'o', 'LineWidth', 2)
%         hold(ax1, 'on')
%         plot(ax1,x1d_ref(i2),p_ref(n,i2),'o', 'LineWidth', 2)    

        fprintf('n=%i, t=%0.4f\n',n,t)            
        
        % domain plot
        plot(ax1,x1d_ref,p_sim(n_i,:),'-', 'LineWidth', 2)
        hold(ax1, 'on')
        plot(ax1,x1d_ref,p_ref(n_i,:),'--', 'LineWidth', 2)
        xline(ax1,r0_pos, '--', 'LineWidth', 2)        
        hold(ax1, 'off')
        xlim(ax1,[x1d_ref(1),x1d_ref(end)])
        ylim(ax1,[-1,1])        
        legend(ax1,'FM-SEM', 'reference', 'receiver')
        xlabel(ax1,'x [m]')
        ylabel(ax1,'Pressure [Pa]', 'Interpreter','latex')
        if write_gif
            title(ax1,sprintf('FM%i-SEM%i \n t=%0.4f sec',data.ppw1,data.ppw2,t))
        else
            title(ax1,sprintf('FM%i-SEM%i',data.ppw1,data.ppw2))
        end
        set(ax1,'FontSize',16)

        % receiver plot
        plot(ax2,tvec(1:n),p_sim(1:n,rec_i),'-', 'LineWidth', 2)
        hold(ax2, 'on')
        plot(ax2,tvec(1:n),p_ref(1:n,rec_i),'--', 'LineWidth', 2)
        hold(ax2, 'off')
        xlim(ax2,[0,iter*s1.domain.dt])
        ylim(ax2,[-0.55,0.55])
        xlabel(ax2,'time [s]')
        ylabel(ax2,sprintf('Pressure [Pa] at x=%0.1f [m]', r0_pos), 'Interpreter','latex')
        legend(ax2,'FM-SEM', 'reference')
        set(ax2,'FontSize',16)
        
        % receiver err plot
        plot(ax3,tvec(1:n),abs(p_sim(1:n,rec_i)-p_ref(1:n,rec_i)),'-', 'LineWidth', 2)
        xlim(ax3,[0,iter*s1.domain.dt])
        ylim(ax3,[0,0.09])
        xlabel(ax3,'time [s]')
        ylabel(ax3,'$\vert p_{FM-SEM} - p_{ref}\vert$ [Pa]', 'Interpreter','latex')
        %ylabel(ax3,'Average droplet velocity [$\mu$/ms]', 'Interpreter','latex')        
        set(ax3,'FontSize',16)

        drawnow
    end

    if mod(n-1,10) == 0 || n==iter
        if write_gif
            writeGif(h, sprintf('%s/coupled_%s_ppw%i_%s_ppw%i.gif',base_path,s1.solver_type,data.ppw1,s2.solver_type,data.ppw2), n==1);
        end
    end
end

if ~write_gif
    %% plot frequency spectrum
    assert(abs(1/s1.domain.dt - 1/s_ref.domain.dt) < eps('double'))
    fs = 1/s1.domain.dt;

    [source, src_x, src_t] = sources.gaussianSourceDiff1D(c, fmax);

    ir = itaAudio();
    ir.samplingRate = fs;
    ir.channelNames{1} = 'Mono channel';
    ir.time = p_sim(:,rec_i);
    ir.trackLength = tmax; %length(ir.time)/fs;

    ir_ref = itaAudio();
    ir_ref.samplingRate = fs;
    ir_ref.channelNames{1} = 'Mono channel';
    ir_ref.time = p_ref(:,rec_i);
    ir_ref.trackLength = tmax;

    src = itaAudio();
    src.samplingRate = fs;
    src.channelNames{1} = 'Mono channel';
    src.time = src_t((1:iter)*s_ref.domain.dt)';
    src.trackLength = tmax;

    hsim = ita_divide_spk(ir, src, 'regularization', [20,fmax]);
    href = ita_divide_spk(ir_ref, src,'regularization', [20,fmax]);

    from = 20;
    plotFrequencyDomain(ax4,hsim.time,fs,fmax,from,"", '-');
    hold(ax4, 'on')
    plotFrequencyDomain(ax4,href.time,fs,fmax,from,"", '--');
    %plotFrequencyDomain(ax4,src.time,fs,fmax,from,"", '--');
    hold(ax4, 'off')
    set(ax4,'FontSize',16)
    ylim(ax4,[-90 -40])
    legend(ax4,'FM-SEM', 'reference', 'location', 'southwest')

    saveas(h,sprintf('%s/FOURIER_SEM_ppw%i_%i_fmax_%i',plots_path,data.ppw1,data.ppw2,fmax),'epsc')
end

%% calculate pulse width
[~,x0_true] = getClosestIndex(x1d_ref, 5.25);
% figure(3)
% ax1 = gca;
% plot(ax1,x1d_ref,p_ref(n,:),'--', 'LineWidth', 2)
% hold(ax1, 'on')
% plot(ax1,x1d_ref,src_x(x1d_ref,0,x0_true)/max(src_x(x1d_ref,0,s_ref.src.x0))/2, '-', 'LineWidth', 2)
% hold(ax1, 'off')
% xlim(ax1,[x1d_ref(1),x1d_ref(end)])
% ylim(ax1,[-1,1])        
% legend(ax1,'reference', 'true source')
% xlabel(ax1,'x [meter]')
% ylabel(ax1,'Pressure [Pa]')
% title(ax1,sprintf('FM %i ppw <-> SEM %i ppw',data.ppw1,data.ppw2))
% set(ax1,'FontSize',16)

[~,indx1_max] = max(src_x(x1d_ref,0,x0_true));
[~,indx1_min] = min(src_x(x1d_ref,0,x0_true));

[~,indx1ref_max] = max(p_ref(n,1:round(length(p_ref(n,:))*2/3)));
[~,indx1ref_min] = min(p_ref(n,1:round(length(p_ref(n,:))*2/3)));

fprintf('Width true: %0.3f\n', abs(x1d_ref(indx1_max) - x1d_ref(indx1_min)))
fprintf('Width ref: %0.3f\n', abs(x1d_ref(indx1ref_max) - x1d_ref(indx1ref_min)))

%% FREQ ERRORS
% calculate errors in frequncy domain
% indxFrom = 1+floor(from/freqs(end) * length(freqs));
% indxTo = min(floor(fn/freqs(end) * length(freqs)),1024/2);
% fROI = freqs(indxFrom:indxTo);
% XXroi_sim = XX_sim(indxFrom:indxTo);
% XXroi_ref = XX_ref(indxFrom:indxTo);
% 
% rms_freq_err = rms(abs(abs(XXroi_sim) - abs(XXroi_ref)));
% rel_freq_err = mean(abs(abs(XXroi_sim) - abs(XXroi_ref))); %./abs(XXroi_ref));
% max_freq_err = max(abs(abs(XXroi_sim) - abs(XXroi_ref)));
% 
% disp('Frequency domaun error measures:')
% fprintf('rms freq err: %0.4f\n', rms_freq_err)
% fprintf('relative freq err: %0.4f dB \n',20*log10(rel_freq_err))
% fprintf('max freq err: %0.4f\n',max_freq_err)
