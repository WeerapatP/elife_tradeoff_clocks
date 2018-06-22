function mutual_information = simulate_millar(signal_min, signal_max, v1, Omega, has_extnoise, has_intnoise)
    % Simulate the Arabidopsis circadian clock by Millar
    % modelling genetic networks with noisy and varied experimental data: the circadian clock in arabidopsis thaliana
    %% to be adjusted
    show_video   = false;
    t_max = 400;
    start_observing = 300;
    hrday = 24.0;
    ext_noise_lvl = 0.5;
    noise_correlation_time = 3/24*hrday;
    
    %% usually fixed
%     dt    = 0.002;
    dt    = 0.01;
    MI_dt = 1;
    vid_dt= MI_dt;
    histogramDx = 0.05;
    times = dt:dt:t_max;
    ensemblesize = 5000;
    minX = inf; minY = inf;
    maxX = -inf; maxY = -inf;
    startrecord = 0;
    ML  = normrnd(0, 0.1, 1, ensemblesize);
    PL  = normrnd(0, 0.1, 1, ensemblesize);
    MT  = normrnd(0, 0.1, 1, ensemblesize);
    PT  = normrnd(0, 0.1, 1, ensemblesize);
    
    video_frame_number = 0;
    HxGivenT=0;
    
    sample_traj_xx = zeros(1, length(times));
    sample_traj_yy = zeros(1, length(times));
    sample_signal  = zeros(1, length(times));
    
    LAmp   = zeros(1, ensemblesize);
    Lones  = ones(1, ensemblesize);
    LNight = zeros(1, ensemblesize);
    
    for i = 1:length(times)
        t = times(i);
        % print t to help me be patient
        if mod(i,  round(100/dt)) == 0
            t
        end
        
        %% light intensity
        if has_extnoise
            change = (rand(1, ensemblesize) < dt/noise_correlation_time);
            LAmp   = LAmp.*(change == 0) + ((1 - ext_noise_lvl) + ext_noise_lvl*rand(1, ensemblesize)).*change;
        else
            LAmp   = ones(1, ensemblesize);
        end
        
        % day or night?
        if (mod(t, hrday) < hrday/2)
            signal = signal_min + (signal_max - signal_min)*LAmp;
        else
            signal = signal_min*ones(1, ensemblesize);
        end
        %% update M, PC, PN
% v1=0.3;
         a=2; g1=0.5; m1=0.4; k1=1;
        p1=0.5; m2=0.6; k2=0.5;
        v2=0.6; b=2; g2=0.1; m3=0.6; k3=1;
        p2=0.3; m4=0.3; k4=1;
        
        PTa      = PT.^a;
        hill_PT   = v1*PTa./(g1^a + PTa);
        ML_decay  = m1*ML./(k1 + ML);
        
        ML_to_PL  = p1*ML;
        PL_decay  = m2*PL./(k2 + PL);
        
        hill_PL   = v2*g2^b./(g2^b + PL.^b);
        MT_decay  = m3*MT./(k3 + MT);
        
        MT_to_PT  = p2*MT;
        PT_decay  = m4*PT./(k4 + PT);
        
        dML      = dt*(signal + hill_PT - ML_decay);
        dPL      = dt*(ML_to_PL - PL_decay);
        dMT      = dt*(hill_PL - MT_decay);
        dPT      = dt*(MT_to_PT - PT_decay);
        if has_intnoise
            dML  = dML + sqrt(abs(dt*hill_PT)/Omega).*normrnd(0, 1, 1, ensemblesize)...
                       + sqrt(abs(dt*ML_decay)/Omega).*normrnd(0, 1, 1, ensemblesize);
            dPL  = dPL + sqrt(abs(dt*ML_to_PL)/Omega).*normrnd(0, 1, 1, ensemblesize)...
                       + sqrt(abs(dt*PL_decay)/Omega).*normrnd(0, 1, 1, ensemblesize);
            dMT  = dMT + sqrt(abs(dt*hill_PL)/Omega).*normrnd(0, 1, 1, ensemblesize)...
                       + sqrt(abs(dt*MT_decay)/Omega).*normrnd(0, 1, 1, ensemblesize);
            dPT  = dPT + sqrt(abs(dt*MT_to_PT)/Omega).*normrnd(0, 1, 1, ensemblesize)...
                       + sqrt(abs(dt*PT_decay)/Omega).*normrnd(0, 1, 1, ensemblesize);
        end
        
        ML      = poslin(ML  + dML);
        PL      = poslin(PL  + dPL);
        MT      = poslin(MT  + dMT);
        PT      = poslin(PT  + dPT);
        
        sample_traj_xx(i)     = ML(1);
        sample_traj_yy(i)     = PL(1);
        sample_signal(i)      = signal(1);
        
        %% Parts of MI calculation
        
        % find max, min
        if(startrecord==0 && (t >= start_observing - hrday*2))
           if(minX > min(ML))
               minX = min(ML);
           end
           if(maxX < max(ML))
               maxX = max(ML);
           end
           if(minY > min(PL))
               minY = min(PL);
           end
           if(maxY < max(PL))
               maxY = max(PL);
           end
        end
        
        % find the right grid for hist
        if((t >= start_observing - dt) && startrecord==0)
            startrecord=1;
            spaceFactor=1.05;
            minX = minX/spaceFactor; minY = minY/spaceFactor;
            maxX = spaceFactor*maxX; maxY = spaceFactor*maxY;
            maxX = maxX + histogramDx;
            maxY = maxY + histogramDx;
            xBin = minX:histogramDx:maxX;
            yBin = minY:histogramDx:maxY;
            numberOfBinsX = length(xBin);
            numberOfBinsY = length(yBin);
            ctr = {yBin xBin};
            N=hist3([PL', ML'],ctr);
            Pxy = zeros(size(N));
            PxyGivenT = zeros(size(N));
        end
        
        % addup Pxy and PxyGivenT
        if((t > t_max - hrday ) && startrecord==1)
            N         = hist3([PL',ML'],ctr);
            N         = N/ensemblesize;
            PxyGivenT = PxyGivenT + N*dt/MI_dt;
            Pxy       = Pxy + N*dt/hrday;
            
            if (show_video && mod(i, round(vid_dt/dt)) == 0)
                video_frame_number = video_frame_number + 1;
                
                figure(237);
                hold all;
                h1 = subplot(2,1,1);
                cla(h1);
                hold all;
                surf(PxyGivenT, 'EdgeColor', 'None');
                title(['time = ', num2str(mod(t, hrday)), ' hrs']);
                xlabel(['[ minX, maxX ] = [ ', num2str(minX), ', ', num2str(maxX), ' ]']);
                ylabel(['[ minY, maxY ] = [ ', num2str(minY), ', ', num2str(maxY), ' ]']);
                view(2);
                
                subplot(2,1,2);
                scatter(mod(t, hrday), sum(-N(N~=0).*log2(N(N~=0))),'o')
                xlabel('time');
                ylabel('entropy');
                hold off;
                view(2);
                set(gcf, 'Position', [0 0 500 800]);
                set(gca,'fontsize',18);
                MM(video_frame_number) = getframe(gcf);
            end
            
            if mod(i, round(MI_dt/dt)) == 0
                HxGivenT = HxGivenT + (-sum(PxyGivenT(PxyGivenT~=0).*log2(PxyGivenT(PxyGivenT~=0)))*MI_dt/hrday);
%                 'check sum'
%                 ['sum(sum(PxyGivenT)) = ', num2str(sum(sum(PxyGivenT)))]
                PxyGivenT = zeros(size(N));
            end
            
        end
    end
    
%     'check sum'
%     ['sum(sum(Pxy)) = ', num2str(sum(sum(Pxy)))]

    %% MI calculation
    Htraj = -sum(Pxy(Pxy~=0).*log2(Pxy(Pxy~=0)));
    mutual_information = Htraj - HxGivenT;
    
    %% period calculation (only works when there is no noise)
    [pks, locs] = findpeaks(sample_traj_xx, times);
    mean_period = mean(diff(locs))
    
    %% make video
    if show_video
        myVideo = VideoWriter(['vs min = ', num2str(signal_min), ',  vs max = ', num2str(signal_max), ', MI = ', num2str(mutual_information), ' van de pol.avi']);
        myVideo.FrameRate = 5;
        open(myVideo);
        writeVideo(myVideo, MM);
        close(myVideo);
    end
    
    %% plot
    figure(21)
    hold all
    starting_index = start_observing/dt;
    plot(sample_traj_xx(starting_index:end), sample_traj_yy(starting_index:end), 'DisplayName', ['signal = [', num2str([signal_min, signal_max]),'], Period = ', num2str(mean_period)],'LineWidth',3);
    xlabel('x');
    ylabel('y');
    title('Simulation of Millar Model');
    legend('-DynamicLegend')
    
    figure(22)
    hold all
    plot(times, sample_signal, times, sample_traj_xx, 'DisplayName',['signal = [', num2str([signal_min, signal_max]),'], Period = ', num2str(mean_period)])
    legend('-DynamicLegend')
%     ending_pos = [sample_traj_xx(end), sample_traj_yy(end)];
%     num2str(ending_pos)
    
end
