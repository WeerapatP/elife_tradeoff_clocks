function mutual_information = simulate_19_equations(signal_min, signal_max, Omega, has_extnoise, has_intnoise, pointattractor)
    % Simulate the Mammalian clock
    %% to be adjusted
    show_video   = false;
    t_max = 700;
    start_observing = 600;
    hrday = 21.5;
    ext_noise_lvl = 0.5;
    noise_correlation_time = 3/24*hrday;
    
    %% usually fixed
%     dt    = 0.002;
    dt    = 0.01;
    MI_dt = 0.5;
    vid_dt= MI_dt;
    histogramDx = 0.01;
    times = dt:dt:t_max;
    ensemblesize = 10000;
    minX = inf; minY = inf;
    maxX = -inf; maxY = -inf;
    startrecord = 0;
    M_P = normrnd(0, 0.1, 1, ensemblesize);
    M_C = normrnd(0, 0.1, 1, ensemblesize);
    M_B = normrnd(0, 0.1, 1, ensemblesize);
    P_C = normrnd(0, 0.1, 1, ensemblesize);
    C_C = normrnd(0, 0.1, 1, ensemblesize);
    P_CP = normrnd(0, 0.1, 1, ensemblesize);
    C_CP = normrnd(0, 0.1, 1, ensemblesize);
    PC_C = normrnd(0, 0.1, 1, ensemblesize);
    PC_N = normrnd(0, 0.1, 1, ensemblesize);
    PC_CP = normrnd(0, 0.1, 1, ensemblesize);
    PC_NP = normrnd(0, 0.1, 1, ensemblesize);
    B_C = normrnd(0, 0.1, 1, ensemblesize);
    B_CP = normrnd(0, 0.1, 1, ensemblesize);
    B_N = normrnd(0, 0.1, 1, ensemblesize);
    B_NP = normrnd(0, 0.1, 1, ensemblesize);
    I_N = normrnd(0, 0.1, 1, ensemblesize);

    
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
        % miu_sP = 1.5;
        miu_sP = signal;
        
        k_1 = 0.4;
        k_2 = 0.2;
        k_3 = 0.4;
        k_4 = 0.2;
        k_5 = 0.4;
        k_6 = 0.2;
        k_7 = 0.5;
        k_8 = 0.1;
        K_AP = 0.7;
%         K_AC = 0.6;
        K_AC = 0.4;
        K_IB = 2.2;
        k_dmb = 0.01;
        k_dmc = 0.01;
        k_dmp = 0.01;
        k_dn = 0.01;
        k_dnc = 0.12;
        K_d = 0.3;
        K_dp = 0.1;
        K_p = 0.1;
        K_mB = 0.4;
        K_mC = 0.4;
        K_mP = 0.31;
        k_sB = 0.12;
        k_sC = 1.6;
        k_sP = 0.6;
        m = 2;
        n = 4;
        V_1B = 0.5; 
        V_1C = 0.6;
        V_1P = 0.4;
        V_1PC = 0.4;
        V_2B = 0.1;
        V_2C = 0.1;
        V_2P = 0.3;
        V_2PC = 0.1;
        V_3B = 0.5;
        V_3PC = 0.4;
        V_4B = 0.2;
        V_4PC = 0.1;
        V_phos = 0.4; %%%%%%this is never used...
        miu_dBC = 0.5;
        miu_dBN = 0.6;
        miu_dCC = 0.7;
        miu_dIN = 0.8;
        miu_dPC = 0.7;
        miu_dPCC = 0.7;
        miu_dPCN = 0.7;
%         miu_mB = 0.8;
        miu_mB = 0.9;
        miu_mC = 1.0;
        miu_mP = 1.1;
        miu_sB = 1.0;
        miu_sC = 1.1;
        if ~pointattractor
            miu_mB = 0.8;
            K_AC = 0.6;
        end
            



        %1
        dM_P = dt*(miu_sP .* (B_N .^n)./ (K_AP ^ n +B_N .^ n) - miu_mP * (M_P)./ (K_mP + M_P) - k_dmp * M_P);
        %2
        dM_C = dt*(miu_sC * (B_N .^n)./ (K_AC ^ n +B_N .^ n) - miu_mC * (M_C)./ (K_mC + M_C) - k_dmc * M_C);
        %3
        dM_B = dt*(miu_sB * (K_IB ^m)./ (K_IB ^ m +B_N .^ m) - miu_mB * (M_B)./ (K_mB + M_B) - k_dmb * M_B);

        %4
        dP_C = dt*(k_sP * M_P - V_1P * P_C ./ (K_p + P_C) + V_2P * P_CP ./ (K_dp +P_CP) + k_4 * PC_C - k_3 * P_C .* C_C - k_dn * P_C);
        %5
        dC_C = dt*(k_sC * M_C - V_1C * C_C ./ (K_p + C_C) + V_2C * C_CP ./ (K_dp +C_CP) + k_4 * PC_C - k_3 * P_C .* C_C - k_dnc * C_C);
        %6
        dP_CP = dt*(V_1P * P_C ./ (K_p + P_C) - V_2P * P_CP ./ (K_dp + P_CP) - miu_dPC * P_CP ./ (K_d + P_CP) - k_dn * P_CP);
        %7
        dC_CP = dt*(V_1C * C_C ./ (K_p + C_C) - V_2C * C_CP ./ (K_dp + C_CP) - miu_dCC * C_CP ./ (K_d + C_CP) - k_dn * C_CP);

        %8
        dPC_C = dt*(-V_1PC * PC_C ./ (K_p + PC_C) + V_2PC * PC_CP ./(K_dp + PC_CP) - k_4 * PC_C + k_3 * P_C .* C_C + k_2 * PC_N - k_1 * PC_C - k_dn * PC_C);
        %9
        dPC_N = dt*(-V_3PC * PC_N ./ (K_p + PC_N) + V_4PC * PC_NP ./(K_dp + PC_NP) - k_2 * PC_N + k_1 * PC_C - k_7 * B_N .* PC_N + k_8 * I_N - k_dn * PC_N);
        %10
        dPC_CP = dt*(V_1PC * PC_C ./ (K_p + PC_C) - V_2PC * PC_CP ./(K_dp + PC_CP) - miu_dPCC * PC_CP ./ (K_d + PC_CP) -k_dn * PC_CP);
        %11
        dPC_NP = dt*(V_3PC * PC_N ./ (K_p + PC_N) - V_4PC * PC_NP ./(K_dp + PC_NP) - miu_dPCN * PC_NP ./ (K_d + PC_NP) -k_dn * PC_NP);

        %12
        dB_C = dt*(k_sB * M_B - V_1B * B_C ./ (K_p + B_C) + V_2B * B_CP ./(K_dp + B_CP) - k_5 * B_C + k_6 * B_N - k_dn * B_C);
        %13
        dB_CP = dt*( V_1B * B_C ./ (K_p + B_C) - V_2B * B_CP ./ (K_dp + B_CP) - miu_dBC * B_CP ./ (K_d + B_CP) - k_dn * B_CP);
        %14
        dB_N = dt*(- V_3B * B_N ./ (K_p + B_N) + V_4B * B_NP ./ (K_dp + B_NP) + k_5 * B_C - k_6 * B_N - k_7 * B_N .* PC_N + k_8 * I_N - k_dn * B_N);
        %15
        dB_NP = dt*(V_3B * B_N ./ (K_p + B_N) - V_4B * B_NP ./ (K_dp + B_NP) - miu_dBN * B_NP ./ (K_d + B_NP) - k_dn * B_NP);
        %16
        dI_N = dt*(- k_8 * I_N + k_7 * B_N .* PC_N - miu_dIN * I_N ./ (K_d + I_N) - k_dn * I_N);


        if has_intnoise
            % add the noise
            dM_P = dM_P +   sqrt(abs(dM_P)/Omega) .*normrnd(0, 1, 1, ensemblesize);
            dM_C = dM_C +   sqrt(abs(dM_C)/Omega)   .*normrnd(0, 1, 1, ensemblesize);
            dM_B = dM_B +   sqrt(abs(dM_C)/Omega)   .*normrnd(0, 1, 1, ensemblesize);
            dP_C = dP_C +    sqrt(abs(dP_C)/Omega) .*normrnd(0, 1, 1, ensemblesize);
            dC_C = dC_C +    sqrt(abs(dC_C)/Omega) .*normrnd(0, 1, 1, ensemblesize);
            dP_CP = dP_CP +   sqrt(abs(dP_CP)/Omega)  .*normrnd(0, 1, 1, ensemblesize);
            dC_CP = dC_CP +    sqrt(abs(dC_CP)/Omega) .*normrnd(0, 1, 1, ensemblesize);
            dPC_C = dPC_C +    sqrt(abs(dPC_C)/Omega) .*normrnd(0, 1, 1, ensemblesize);
            dPC_N = dPC_N +    sqrt(abs(dPC_N)/Omega) .*normrnd(0, 1, 1, ensemblesize);
            dPC_CP = dPC_CP +    sqrt(abs(dPC_CP)/Omega) .*normrnd(0, 1, 1, ensemblesize);
            dPC_NP = dPC_NP +    sqrt(abs(dPC_NP)/Omega) .*normrnd(0, 1, 1, ensemblesize);
            dB_C = dB_C +   sqrt(abs(dB_C)/Omega)  .*normrnd(0, 1, 1, ensemblesize);
            dB_CP = dB_CP +    sqrt(abs(dB_CP)/Omega) .*normrnd(0, 1, 1, ensemblesize);
            dB_N = dB_N +   sqrt(abs(dB_N)/Omega)  .*normrnd(0, 1, 1, ensemblesize);
            dB_NP = dB_NP +   sqrt(abs(dB_NP)/Omega)  .*normrnd(0, 1, 1, ensemblesize);
            dI_N = dI_N +    sqrt(abs(dI_N)/Omega)  .*normrnd(0, 1, 1, ensemblesize);
        end


        % updating the values
        M_P =  poslin(dM_P + M_P);
        M_C =  poslin(dM_C +  M_C);
        M_B =  poslin(dM_B + M_B);
        P_C =  poslin(dP_C +  P_C);
        C_C =  poslin(dC_C + C_C);
        P_CP =  poslin(dP_CP + P_CP);
        C_CP =  poslin(dC_CP + C_CP);
        PC_C =  poslin(dPC_C + PC_C);
        PC_N =  poslin(dPC_N + PC_N);
        PC_CP =  poslin(dPC_CP + PC_CP);
        PC_NP =  poslin(dPC_NP + PC_NP);
        B_C =  poslin(dB_C + B_C);
        B_CP =  poslin(dB_CP + B_CP);
        B_N =  poslin(dB_N + B_N);
        B_NP =  poslin(dB_NP + B_NP);
        I_N =  poslin(dI_N +I_N);

        xx = M_B;
        yy = M_P;
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
