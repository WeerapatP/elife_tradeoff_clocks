function mutual_information = simulate_repressilator(signal_min, signal_max, alpha, n, Omega, has_extnoise, has_intnoise)
    % Simulate the repressilator equations 
    % Bu?e, O., Kuznetsov, A., and Pérez, R. A. (2009). EXISTENCE OF LIMIT CYCLES IN THE REPRESSILATOR EQUATIONS
    %% to be adjusted
    show_video   = false;
    t_max = 500;
    start_observing = 400;
    hrday = 4;
    ext_noise_lvl = 0.5;
    noise_correlation_time = 3/24*hrday;
    
    %% usually fixed
    dt    = 0.01; 
    MI_dt = 0.2;
    vid_dt= MI_dt;
    histogramDx = 0.2;
    times = dt:dt:t_max;
    ensemblesize = 2000;
    minX = inf; minY = inf;
    maxX = -inf; maxY = -inf;
    startrecord = 0;
    xx  = normrnd(0, 1, 1, ensemblesize);
    yy  = normrnd(0, 1, 1, ensemblesize);
    zz  = normrnd(0, 1, 1, ensemblesize);
    
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
%         alpha    = 10;
        hill_y   = alpha*(signal./(1+yy.^n));
        hill_z   = alpha./(1+zz.^n);
        hill_x   = alpha./(1+xx.^n);
        dxx      = dt*(hill_y - xx);
        dyy      = dt*(hill_z - yy);
        dzz      = dt*(hill_x - zz);
        if has_intnoise
            dxx  = dxx + sqrt(abs(dt*hill_y)/Omega).*normrnd(0, 1, 1, ensemblesize)...
                       + sqrt(abs(dt*xx)/Omega).*normrnd(0, 1, 1, ensemblesize);
                   
            dyy  = dyy + sqrt(abs(dt*hill_z)/Omega).*normrnd(0, 1, 1, ensemblesize)...
                       + sqrt(abs(dt*yy)/Omega).*normrnd(0, 1, 1, ensemblesize);
                   
            dzz  = dzz + sqrt(abs(dt*hill_x)/Omega).*normrnd(0, 1, 1, ensemblesize)...
                       + sqrt(abs(dt*zz)/Omega).*normrnd(0, 1, 1, ensemblesize);
        end
        
        xx      = poslin(xx  + dxx);
        yy      = poslin(yy  + dyy);
        zz      = poslin(zz  + dzz);
        
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
    title('Simulation of Goldbeter Model');
    legend('-DynamicLegend')
    
    figure(22)
    plot(times, sample_signal, times, sample_traj_xx)
%     ending_pos = [sample_traj_xx(end), sample_traj_yy(end)];
%     num2str(ending_pos)
    
end

%         Rungekutta
%         M1  =  M;
%         PC1 = PC;
%         PN1 = PN;
%         dMdt1   = vs*(KI^n./(KI^n+PN1.^n)) - vm*(M1./(Km + M1));
%         dPCdt1  = - k1*PC1 + k2*PN1 + ks*M1 - vd*(PC1./(Kd + PC1));
%         dPNdt1  = + k1*PC1 - k2*PN1;
% 
%         M2  =  M + dMdt1*dt/2;
%         PC2 = PC + dPCdt1*dt/2;
%         PN2 = PN + dPNdt1*dt/2;
%         dMdt2   = vs*(KI^n./(KI^n+PN2.^n)) - vm*(M2./(Km + M2));
%         dPCdt2  = - k1*PC2 + k2*PN2 + ks*M2 - vd*(PC2./(Kd + PC2));
%         dPNdt2  = + k1*PC2 - k2*PN2;
% 
%         M3  =  M + dMdt2*dt/2;
%         PC3 = PC + dPCdt2*dt/2;
%         PN3 = PN + dPNdt2*dt/2;
%         dMdt3   = vs*(KI^n./(KI^n+PN3.^n)) - vm*(M3./(Km + M3));
%         dPCdt3  = - k1*PC3 + k2*PN3 + ks*M3 - vd*(PC3./(Kd + PC3));
%         dPNdt3  = + k1*PC3 - k2*PN3;
% 
%         M4  =  M + dMdt3*dt;
%         PC4 = PC + dPCdt3*dt;
%         PN4 = PN + dPNdt3*dt;
%         dMdt4   = vs*(KI^n./(KI^n+PN4.^n)) - vm*(M4./(Km + M4));
%         dPCdt4  = - k1*PC4 + k2*PN4 + ks*M4 - vd*(PC4./(Kd + PC4));
%         dPNdt4  = + k1*PC4 - k2*PN4;
% 
%         M  = M  + dt*(dMdt1  + 2*dMdt2  + 2*dMdt3  + dMdt4)/6;
%         PC = PC + dt*(dPCdt1 + 2*dPCdt2 + 2*dPCdt3 + dPCdt4)/6;
%         PN = PN + dt*(dPNdt1 + 2*dPNdt2 + 2*dPNdt3 + dPNdt4)/6;
