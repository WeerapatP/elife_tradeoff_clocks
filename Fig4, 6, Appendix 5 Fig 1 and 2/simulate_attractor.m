function [mid_entropy, max_entropy, mutual_information] = simulate_attractor(alpha, r0, Lmax, ensemblesize, dt, noise_correlation_time, sigma, ext_noise_lvl, finalDay)
    % Simulate Attractor System
    %
    % Weerapat Pittayakanchit
    
    DEBUG = true;
    DEBUG_day = 3;
    showFigure = false;
    showVideo  = false;
    if ext_noise_lvl == 0
        hasNoise = false;
    else
        hasNoise = true; 
    end

    n_days = finalDay;% Total number of days, including before genesis
    hd=12;%12 hour daylight
    frac_day=hd/24;
    
    % Fixed Number
    histogramDx = 0.01;
% 	ctrFixed  = {-2*r0:histogramDx:2*r0 2*(-r0):histogramDx:2*(Lmax+r0)};
    
    omega0=2*pi;  % this is omega of night, changing it doesn't change the fact
    % that the period of a night limit cycle is 24 hours.
    sqrtdt=sqrt(dt);
    MI_dt = 0.02;

    %%%%%%%%set initial ensemble on the ring of darkness...
    angles = linspace(0,2*pi,ensemblesize);
    if alpha > 0
        xx = r0*cos(angles');
        yy = r0*sin(angles');
    else
        xx = zeros(ensemblesize, 1);
        yy = zeros(ensemblesize, 1);        
    end
    
    DayLight = ones(ensemblesize,1);
    LDay = (rand(ensemblesize,1)>0.5);
    LAmp = rand(ensemblesize,1);
    LNight = zeros(ensemblesize,1);
    
        
    if DEBUG
        sample_ent = zeros(length(n_days - DEBUG_day:dt:n_days), 1);
        sample_L   = zeros(length(n_days - DEBUG_day:dt:n_days), 1);
        sample_r   = zeros(length(n_days - DEBUG_day:dt:n_days), 1);
        sample_x   = zeros(length(n_days - DEBUG_day:dt:n_days), 1);
    end
    
    %initial state generated, below starts the day-night cycle
    startrecord=0;
    HxGivenT=0;
    
    counter_for_MI = 1;
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
    
    for time=0:dt:n_days
        
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
        
        if(startrecord==0 && (time > n_days - DEBUG_day - 2))
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
        
        if((time > n_days - DEBUG_day - 1) && startrecord==0)
            startrecord=1;
            spaceFactor=1.05;
            minX = spaceFactor*minX; minY = spaceFactor*minY;
            maxX = spaceFactor*maxX; maxY = spaceFactor*maxY;
            maxX = maxX + histogramDx;
            maxY = maxY + histogramDx;
            xBin = minX:histogramDx:maxX;
            yBin = minY:histogramDx:maxY;
            numberOfBinsX = length(xBin);
            numberOfBinsY = length(yBin);
            ctr = {yBin xBin};
%             ctr = ctrFixed;
            N=hist3([yy,xx],ctr);
            Pxy = zeros(size(N));
            PxyGivenT = zeros(size(N));
        end
        
        if DEBUG && time >= n_days - DEBUG_day
            index_signals = index_signals + 1;

            sample_ent(index_signals) = log2(max(eig(cov([xx yy]))));
            sample_L(index_signals) = L(1);
            sample_r(index_signals) = sqrt((xx(1)-Lmax/2).^2 + yy(1).^2);
            sample_x(index_signals) = xx(1);
        end

        %here we calculate
        if((time > n_days - 1) && startrecord==1)
            
            N=hist3([yy,xx],ctr);
            N=N/ensemblesize;
            PxyGivenT = PxyGivenT + N*dt/MI_dt;
            Pxy=Pxy+N*dt;
            
            counter_for_video = counter_for_video + 1;
            if(showVideo && counter_for_video >= max_counter_video)
                counter_for_video = 0;
                video_frame_number = video_frame_number + 1;
                
                figure(237);
                hold all;
                h1 = subplot(2,1,1);
                cla(h1);
                hold all;
                ttt = linspace(0,2*pi);
                xCenter1 = numberOfBinsX*(0-minX)/(maxX-minX);    yCenter1 = numberOfBinsY/2;
                xCenter2 = numberOfBinsX*(Lmax-minX)/(maxX-minX);    yCenter2 = numberOfBinsY/2;
                radius = r0/histogramDx;
                
                if alpha > 0
                    if rem(time,1) < 0.5
                        plot(1 + xCenter1 + radius*cos(ttt), 1 + yCenter1 + radius*sin(ttt),'Color','red','LineWidth',1)
                        plot(1 + xCenter2 + radius*cos(ttt), 1 + yCenter2 + radius*sin(ttt),'Color','green','LineWidth',1)
                    else
                        plot(1 + xCenter1 + radius*cos(ttt), 1 + yCenter1 + radius*sin(ttt),'Color','green','LineWidth',1)
                        plot(1 + xCenter2 + radius*cos(ttt), 1 + yCenter2 + radius*sin(ttt),'Color','red','LineWidth',1)
                    end
                else
                    radius = 0.1;
                    if rem(time,1) < 0.5
                        plot(1 + xCenter1 + radius*cos(ttt), 1 + yCenter1 + radius*sin(ttt),'Color','red','LineWidth',3)
                        plot(1 + xCenter2 + radius*cos(ttt), 1 + yCenter2 + radius*sin(ttt),'Color','green','LineWidth',3)
                    else
                        plot(1 + xCenter1 + radius*cos(ttt), 1 + yCenter1 + radius*sin(ttt),'Color','green','LineWidth',3)
                        plot(1 + xCenter2 + radius*cos(ttt), 1 + yCenter2 + radius*sin(ttt),'Color','red','LineWidth',3)
                    end
                end
                surf(PxyGivenT, 'EdgeColor', 'None');
                title(['time = ', num2str(24*mod(time,1)), ' hrs, H(X|T) = ', num2str((-sum(PxyGivenT(PxyGivenT~=0).*log2(PxyGivenT(PxyGivenT~=0)))))]);
%                 xlabel(['[ minX, maxX ] = [ ', num2str(minX), ', ', num2str(maxX), ' ]']);
%                 ylabel(['[ minY, maxY ] = [ ', num2str(minY), ', ', num2str(maxY), ' ]']);
                xlabel('x')
                ylabel('y')
                view(2);
                
                subplot(2,1,2);
                scatter(mod(time, 1), sum(-N(N~=0).*log2(N(N~=0))),'o')
                xlabel('time');
                ylabel('entropy');
                hold off;
                view(2);
                set(gcf, 'Position', [0 0 500 800]);
                set(gca,'fontsize',18);
                MM(video_frame_number) = getframe(gcf);
            end
            
            counter_for_MI = counter_for_MI + 1;
            if(counter_for_MI >= max_counter_MI)
                counter_for_MI = 0;
                HxGivenT = HxGivenT + (-sum(PxyGivenT(PxyGivenT~=0).*log2(PxyGivenT(PxyGivenT~=0)))*MI_dt);
                PxyGivenT = zeros(size(N));
            end
            
        end

    end
    
    if showFigure
        figure()
        surf(Pxy, 'EdgeColor', 'none')
        view(2)
        title(['Relaxation Time = ', num2str(20/alpha), ' hrs']);
        xlabel(['[ minX, maxX ] = [ ', num2str(minX), ', ', num2str(maxX), ' ]']);
        ylabel(['[ minY, maxY ] = [ ', num2str(minY), ', ', num2str(maxY), ' ]']);
        ttt = linspace(0,2*pi);
        hold all;
        xCenter1 = numberOfBinsX*(0-minX)/(maxX-minX);    yCenter1 = numberOfBinsY*(maxY-0)/(maxY-minY);
        xCenter2 = numberOfBinsX*(Lmax-minX)/(maxX-minX); yCenter2 = numberOfBinsY*(maxY-0)/(maxY-minY);
        radius = numberOfBins*r0/(maxX-minX);
        plot(1 + xCenter1 + radius*cos(ttt), 1 + yCenter1 + radius*sin(ttt),'Color','red','LineWidth',3)
        plot(1 + xCenter2 + radius*cos(ttt), 1 + yCenter2 + radius*sin(ttt),'Color','red','LineWidth',3)
        hold off;
    end
    
    Htraj = -sum(Pxy(Pxy~=0).*log2(Pxy(Pxy~=0)));
    mutual_information= Htraj - HxGivenT;
    
    if showVideo
        myVideo = VideoWriter(['noise time = ', num2str(noise_correlation_time*24), ' hrs, r0 = ', num2str(r0), ', L = ', num2str(Lmax), ', alpha = ', num2str(alpha), ', hasNoise = ',  num2str(hasNoise), ', MI = ', num2str(mutual_information), '.avi']);
        myVideo.FrameRate = 5;
        open(myVideo);
        writeVideo(myVideo, MM);
        close(myVideo);
    end
   
    
    if DEBUG
        figure()
        plot(n_days - DEBUG_day:dt:n_days, sample_ent)
        xlabel('time')
%         ylabel('Entropy')
        ylabel('Entropy ~ log2(Variance)')
        title(['r0 = ', num2str(r0), ', L = ', num2str(Lmax), ', alpha = ', num2str(alpha), ', ext noise lvl = ', num2str(ext_noise_lvl)])
        set(gca,'fontsize',12);
%         figure(102)
%         plot(24*(0:dt:DEBUG_day), sample_L, '-', 24*(0:dt:DEBUG_day),sample_r-r0, '-')
%         xlabel('Time (hrs)')
%         ylabel('Length')
%         legend('L','dr');
%         title(['r0 = ', num2str(r0), ', L = ', num2str(Lmax), ', alpha = ', num2str(alpha)])
%         set(gca,'fontsize',12);
%         saveas(gcf,['plot r and L vs time for r0 = ', num2str(r0), ', L = ', num2str(Lmax), ', alpha = ', num2str(alpha), '.pdf'],'pdf')
        max_entropy = max(sample_ent);
        min_entropy = min(sample_ent);
        mid_entropy = sample_ent(end);
        
        figure()
        xlabel('time')
        ylabel('x')
        plot(24*(0:dt:DEBUG_day), sample_x)

        sigma_e = max_entropy - min_entropy;
        entropy_drop_dawn = mid_entropy - min_entropy;
        entropy_drop_dusk = max_entropy - mid_entropy;
    end
end
