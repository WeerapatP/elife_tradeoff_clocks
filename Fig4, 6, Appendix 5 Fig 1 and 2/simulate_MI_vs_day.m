function [entropies, mutual_informations] = simulate_MI_vs_day(alpha, r0, Lmax, ensemblesize, dt, noise_correlation_time, sigma, ext_noise_lvl, finalDay)

    % simulate attractor system: modified to show more information than the
    % function simulate_attractor
    %
    % Weerapat Pittayakanchit

    DEBUG = true;
    DEBUG_day = 5;
    showFigure = false;
    showVideo  = false;
    
    if ext_noise_lvl == 0
        hasNoise = false;
    else
        hasNoise = true; 
    end

    n_days=finalDay;% Total number of days, including before genesis
    hd=12;%12 hour daylight
    frac_day=hd/24;
    
    % Fixed Number
    histogramDx = 0.01;
    
    omega0=2*pi;  % this is omega of night, changing it doesn't change the fact
    % that the period of a night limit cycle is 24 hours.
    sqrtdt=sqrt(dt);
    MI_dt = 0.02;

    %%%%%%%%set initial ensemble on the ring of darkness...
    angles = linspace(0,2*pi,ensemblesize);
    if alpha > 0
        xx = r0*cos(angles');
        yy = r0*sin(angles');
        ctrFixed  = {-2*r0:histogramDx:2*r0 2*(-r0):histogramDx:(3*r0+Lmax)};
    else
        xx = zeros(1,length(ensemblesize));
        yy = zeros(1,length(ensemblesize));      
        ctrFixed  = {-2*Lmax:histogramDx:2*Lmax 2*(-Lmax):histogramDx:(3*Lmax)};
    end
    
    DayLight = ones(ensemblesize,1);
    LDay = (rand(ensemblesize,1)>0.5);
    LAmp = rand(ensemblesize,1);
    LNight = zeros(ensemblesize,1);
    
        
    if DEBUG
        sample_ent = zeros(length(n_days - DEBUG_day:dt:n_days), 1);
        sample_L   = zeros(length(n_days - DEBUG_day:dt:n_days), 1);
        sample_r   = zeros(length(n_days - DEBUG_day:dt:n_days), 1);
    end
    
    %initial state generated, below starts the day-night cycle
    startrecord=0;
    max_counter_MI = floor(MI_dt/dt);

    minX = inf; minY = inf;
    maxX = -inf; maxY = -inf;
    
    % set up variables for plotting video
    video_dt          = MI_dt;
    max_counter_video = floor(video_dt/dt);
    counter_for_video = 1;
    video_frame_number = 0;
    
    % sample the signal
    index_signals = 0;
    if showVideo
        figure(237)
        clf;
    end
    
    time_index = 0;
    entropies = zeros(1, length(0:dt:n_days));
    starting_day_dt     = 0.2;
    mutual_informations = zeros(1,n_days);
    
    for time = 0:dt:n_days
        
        % Set 0 at the beginning of the day to find new MI each day
        if mod(time,1) < dt/2
            % Calculate MI at the end of the previous day
            if time > 0
                Htraj = -sum(Pxy(Pxy~=0).*log2(Pxy(Pxy~=0)));
                mutual_informations(round(time)) = Htraj - HxGivenT;
%                 checkSum = sum(sum(Pxy))
            end
            
            % Set 0 after calculating MI
            N=hist3([yy,xx],ctrFixed);
            Pxy = zeros(size(N));
            PxyGivenT = zeros(size(N));
            HxGivenT=0;
            counter_for_MI = max_counter_MI;
        end
        
        % ext_noise_lvl is assummed to be from 0 to 1
        if hasNoise
            change = (rand(ensemblesize,1) < dt/noise_correlation_time);
            LAmp = LAmp.*(change == 0) + ((1 - ext_noise_lvl) + ext_noise_lvl*rand(ensemblesize,1)).*change;
        else
            LAmp = DayLight;
        end
        
        if (rem(time,1) < frac_day)
            L = Lmax*LAmp;
        else
            L = Lmax*LNight;
        end


        dxx = dt*( alpha*(1-sign(alpha)*((xx-L).^2 +yy.^2)/r0^2).*(xx-L) - yy*omega0);
        dyy = dt*( alpha*(1-sign(alpha)*((xx-L).^2 + yy.^2)/r0^2).*yy + (xx-L)*omega0);
        
        xx = xx + dxx;
        yy = yy + dyy;
        
        if sigma ~= 0
            xx=xx + normrnd(0.0,sigma, ensemblesize, 1)*sqrtdt;
            yy=yy + normrnd(0.0,sigma, ensemblesize, 1)*sqrtdt;
        end
        
        % Calculate entropy for each time step
        N=hist3([yy,xx],ctrFixed);
        N=N/ensemblesize;     % probability
        
        time_index = time_index + 1;
        entropies(time_index) = -sum(N(N~=0).*log2(N(N~=0)));
        
        Pxy = Pxy+N*dt;
        PxyGivenT = PxyGivenT + N*dt/MI_dt;

        if(counter_for_MI == max_counter_MI)
            counter_for_MI = 0;
            HxGivenT = HxGivenT + (-sum(PxyGivenT(PxyGivenT~=0).*log2(PxyGivenT(PxyGivenT~=0)))*MI_dt);
            PxyGivenT = zeros(size(N));
        end
        counter_for_MI = counter_for_MI + 1;
    
    end
    
    if showVideo
        myVideo = VideoWriter(['noise time = ', num2str(noise_correlation_time*24), ' hrs, r0 = ', num2str(r0), ', L = ', num2str(Lmax), ', alpha = ', num2str(alpha), ', hasNoise = ',  num2str(hasNoise), ', MI = ', num2str(mutual_informations(end)), '.avi']);
        myVideo.FrameRate = 5;
        open(myVideo);
        writeVideo(myVideo, MM);
        close(myVideo);
    end
end
   
