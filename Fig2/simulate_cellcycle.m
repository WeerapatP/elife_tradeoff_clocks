function mutual_information = simulate_cellcycle(signal_min, signal_max, Omega, has_extnoise, has_intnoise)
    % Simulate cdc2-cyclin cell cycle by Goldbeter
    %% to be adjusted
    show_video   = false;
    t_max = 1000;
    start_observing = 900;
    hrday = 35;
    ext_noise_lvl = 0.7;
    noise_correlation_time = 3/24*hrday;
    
    %% usually fixed
%     dt    = 0.002;
    dt    = 0.1;
    MI_dt = 0.5;
    vid_dt= MI_dt;
    histogramDx = 0.005;
    times = dt:dt:t_max;
    ensemblesize = 12000;
    minX = inf; minY = inf;
    maxX = -inf; maxY = -inf;
    startrecord = 0;
    
    C = normrnd(0, 0, 1, ensemblesize)+0.01;
    M = normrnd(0, 0, 1, ensemblesize)+0.01;
    X = normrnd(0, 0, 1, ensemblesize)+0.01;

    
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
        % print t to h23elp me be patient
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
        %miu_i =0.075 ;%original paper use 0.025; range for point attractor 0.001-0.0105
%         miu_d = 0.25;
%         K_d = 0.02;
%         k_d = 0.01;
%         V_M1 = 3;
%         V_2 = 1.5;
%         V_M3 = 1;
%         V_4 = 0.5;
%         K_C = 0.5;
%         K_1 = 0.005;
%         K_2 = 0.005;
%         K_3 = 0.005;
%         K_4 = 0.005;

        K_1 = 0.1;
        K_2 = 0.1;
        K_3 = 0.1;
        K_4 = 0.1;
        K_C = 0.3;
        V_M1= 0.5;
        V_2 = 0.167;
        V_M3= 0.2;
        V_4 = 0.1;
        miu_d = 0.1;
        K_d = 0.02;
        k_d = 3.33/1000;
        miu_i = signal;
        


        V_1 = C .* V_M1 ./ (K_C * Omega + C);
        V_3 = M .* V_M3;
        CBirth = miu_i * Omega;
        CDecay = k_d * C;
        CDecay_X = miu_d *Omega * X .* C ./(K_d*Omega + C);
        MBirth = V_1 .* (1 - M) ./ (K_1 + 1 - M);
        MDecay = V_2 * M ./ (K_2 + M);
        XBirth = V_3 .* (1-X) ./ (K_3 + 1 - X);
        XDecay = V_4 * X ./ (K_4 + X);
        
        dC = dt* (CBirth - CDecay_X - CDecay);
        dM = dt* (MBirth - MDecay);
        dX = dt* (XBirth - XDecay);

        % add the noise
        if (has_intnoise)
        dC = dC +   sqrt(dt*abs(CBirth))   .*normrnd(0, 1, 1, ensemblesize)...
                +   sqrt(dt*abs(CDecay_X)) .*normrnd(0, 1, 1, ensemblesize)...
                +   sqrt(dt*abs(CDecay))   .*normrnd(0, 1, 1, ensemblesize);
        dM = dM +   sqrt(dt*abs(MBirth)/Omega)   .*normrnd(0, 1, 1, ensemblesize)...
                +   sqrt(dt*abs(MDecay)/Omega)   .*normrnd(0, 1, 1, ensemblesize);
        dX = dX +   sqrt(dt*abs(XBirth)/Omega)   .*normrnd(0, 1, 1, ensemblesize)...
                +   sqrt(dt*abs(XDecay)/Omega)   .*normrnd(0, 1, 1, ensemblesize);
        end



        % updating the values
        C = poslin(C + dC);
        M = poslin(M + dM);
        X = poslin(X + dX);
        M(M>1) = 1;
        X(X>1) = 1;
 

        xx = M;
        yy = X;
        sample_traj_xx(i)     = xx(1);
        sample_traj_yy(i)     = yy(1);
        sample_signal(i)      = signal(1);
        
        %% Parts of MI calculation
        
        % find max, min
        if(startrecord==0 && (t >= start_observing - hrday*2))
           if(minX > min(xx))
               minX = min(xx);
           end
           if(maxX < max(xx))
               maxX = max(xx);
           end
           if(minY > min(yy))
               minY = min(yy);
           end
           if(maxY < max(yy))
               maxY = max(yy);
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
            N=hist3([yy', xx'],ctr);
            Pxy = zeros(size(N));
            PxyGivenT = zeros(size(N));
        end
        
        % addup Pxy and PxyGivenT
        if((t > t_max - hrday ) && startrecord==1)
            N         = hist3([yy',xx'],ctr);
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
    
    starting_index = start_observing/dt;
    [pks, locs] = findpeaks(sample_traj_xx(starting_index:end), times(starting_index:end));
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
