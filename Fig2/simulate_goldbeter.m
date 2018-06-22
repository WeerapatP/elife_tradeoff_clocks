function mutual_information = simulate_goldbeter(vs_min, vs_max, Omega, has_extnoise, has_intnoise)
    % Simulate Neurospora and Drosophilla circadian clocks by Goldbeter
    show_video   = true;
    
    KI    = 1;
    n     = 4;
    vm    = 0.505;
    Km    = 0.5;
    ks    = 0.5;
    vd    = 1.4;
    Kd    = 0.13;
    k1    = 0.1;
%     k1    = 0.5;
    k2    = 0.6;
    dt    = 0.01; % in hour
    MI_dt = 1;
    vid_dt= 1;
    t_max = 1000;
    start_observing = 900;
%     t_max = 2000;
%     start_observing = 1900;
    
    hrday = 18;
    noise_correlation_time = 3;
    histogramDx = 0.05;
    ext_noise_lvl = 0.3;
    times = dt:dt:t_max;
    ensemblesize = 10000;
    minX = inf; minY = inf;
    maxX = -inf; maxY = -inf;
    startrecord = 0;
%     M  = zeros(1, ensemblesize);
%     PC = zeros(1, ensemblesize);
%     PN = zeros(1, ensemblesize);
    M  = Omega*abs(normrnd(0, 1, 1, ensemblesize));
    PC = Omega*abs(normrnd(0, 1, 1, ensemblesize));
    PN = Omega*abs(normrnd(0, 1, 1, ensemblesize));
    
    video_frame_number = 0;
    HxGivenT=0;
    
    sample_traj_M  = zeros(1, length(times));
    sample_traj_PC = zeros(1, length(times));
    sample_traj_PN = zeros(1, length(times));
    sample_vs      = zeros(1, length(times));
    
    LAmp   = zeros(1, ensemblesize);
    Lones  = ones(1, ensemblesize);
    LNight = zeros(1, ensemblesize);
    
    for i = 1:length(times)
        t = times(i);
%         if mod(i,  round(100/dt)) == 0
%             t
%         end
        
        %% light intensity
        if has_extnoise
            change = (rand(1, ensemblesize) < dt/noise_correlation_time);
            LAmp   = LAmp.*(change == 0) + ((1 - ext_noise_lvl) + ext_noise_lvl*rand(1, ensemblesize)).*change;
        else
            LAmp   = ones(1, ensemblesize);
        end
        
        if (mod(t, hrday) < hrday/2)
            vs = vs_min + (vs_max - vs_min)*LAmp;
        else
            vs = vs_min*ones(1, ensemblesize);
        end
        %% update M, PC, PN
        Mbirth  = (vs*Omega).*((KI*Omega)^n./((KI*Omega)^n+PN.^n));
        Mdecay  = (vm*Omega)*(M./((Km*Omega) + M));
        PCbirth = ks*M;
        PCdecay = (vd*Omega)*(PC./((Kd*Omega) + PC));
        PCtoPN  = k1*PC;
        PNtoPC  = k2*PN;
        
        if has_intnoise
            dM      = dt*(Mbirth - Mdecay)...
                        + sqrt(dt*Mbirth).* normrnd(0, 1, 1, ensemblesize)...
                        + sqrt(dt*Mdecay).* normrnd(0, 1, 1, ensemblesize);
            dPC     = dt*(  PNtoPC - PCtoPN + PCbirth - PCdecay)...
                        + sqrt(dt*PNtoPC).* normrnd(0, 1, 1, ensemblesize)...
                        + sqrt(dt*PCtoPN).* normrnd(0, 1, 1, ensemblesize)...
                        + sqrt(dt*PCbirth).*normrnd(0, 1, 1, ensemblesize)...
                        + sqrt(dt*PCdecay).*normrnd(0, 1, 1, ensemblesize);
            dPN     = dt*(- PNtoPC + PCtoPN)...
                        + sqrt(dt*PNtoPC).* normrnd(0, 1, 1, ensemblesize)...
                        + sqrt(dt*PCtoPN).* normrnd(0, 1, 1, ensemblesize);
        else
            dM      = dt*(  Mbirth - Mdecay);
            dPC     = dt*(  PNtoPC - PCtoPN + PCbirth - PCdecay);
            dPN     = dt*(- PNtoPC + PCtoPN);
        end
        
        M      = poslin(M  + dM);
        PC     = poslin(PC + dPC);
        PN     = poslin(PN + dPN);
        sample_traj_M(i)  =  M(1);
        sample_traj_PC(i) = PC(1);
        sample_traj_PN(i) = PN(1);
        sample_vs(i)      = vs(1);
        
        %% Parts of MI calculation
        xx = M/Omega;
        yy = PC/Omega;
        
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
    [pks, locs] = findpeaks(sample_traj_M, times);
    mean_period = mean(diff(locs));
    
    %% make video
    if show_video
        myVideo = VideoWriter(['vs min = ', num2str(vs_min), ',  vs max = ', num2str(vs_max), ', Omega = ', num2str(Omega), ', MI = ', num2str(mutual_information), '.avi']);
        myVideo.FrameRate = 5;
        open(myVideo);
        writeVideo(myVideo, MM);
        close(myVideo);
    end
    
    %% plot
    figure(21)
    hold all
    starting_index = start_observing/dt;
    plot(sample_traj_M(starting_index:end)/Omega, sample_traj_PC(starting_index:end)/Omega, 'DisplayName', ['vs = [', num2str([vs_min, vs_max]),'], Period = ', num2str(mean_period)],'LineWidth',3);
    xlabel('M concentration');
    ylabel('PC concentration');
    title('Simulation of Goldbeter Model');
    
    figure(22)
    plot(times, sample_traj_M)
    
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
