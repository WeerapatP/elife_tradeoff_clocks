alpha = 5;
gamma = 0;
r0    = 0.5;
ensemblesize = 10000;
dt    = 0.001;
noise_correlation_time = 0.1;
finalDay = 150;
numberOfData = 12;

alphas          = alpha*ones(1, numberOfData);
alphas(end)     = -alpha;

r0s             = r0*ones(1, numberOfData);
r0s(end)        = 1000;

L_attractor     = 1;
Lmaxs           = logspace(log10(0.03),log10(1), numberOfData - 1);
Lmaxs           = [Lmaxs L_attractor];

% use this for external noise Fig. 6b
% sigmas          = linspace(0,0, numberOfData);
% ext_noise_lvls  = linspace(1,1, numberOfData);

% use this for internal noise Fig. 6c
sigmas          = linspace(0.1,0.1, numberOfData);
ext_noise_lvls  = linspace(0,0, numberOfData);

mutual_informations = zeros(length(sigmas), length(1:finalDay));
entropies           = zeros(length(sigmas), length(0:dt:finalDay));

% ext_noise_lvls & sigmas
parfor i = 1:numberOfData
    [entropies(i,:) , mutual_informations(i,:)] = simulate_MI_vs_day(alpha, r0, Lmaxs(i), ensemblesize, dt, noise_correlation_time, sigmas(i), ext_noise_lvls(i), finalDay);
end

figure()
for i = 1:numberOfData
    plot(1:finalDay, mutual_informations(i, :), 'DisplayName', sprintf('sigma = %0.2f, ext noise lvls = %0.2f, L = %0.2f', sigmas(i), ext_noise_lvls(i), Lmaxs(i)))
    xlabel('# day')
    ylabel('MI')
    hold all;
end
hold off;
title(sprintf('r0 = %0.2f, alpha = %0.2f', r0, alpha));
legend('-DynamicLegend')

figure()
for i = 1:length(sigmas)
    plot(1:finalDay-1, diff(mutual_informations(i, :)), 'DisplayName', sprintf('sigma = %0.2f, ext noise lvls = %0.2f, L = %0.2f', sigmas(i), ext_noise_lvls(i), Lmaxs(i)))
    xlabel('# day')
    ylabel('MI change')
%     title(sprintf('r0 = %0.2f, L = %0.2f, alpha = %0.2f, sigma = %0.2f, ext noise lvl = %0.2f', r0, Lmax, alpha, sigma, ext_noise_lvl))
    hold all;
end
hold off;
title(sprintf('r0 = %0.2f, alpha = %0.2f', r0, alpha));
legend('-DynamicLegend')

time_almost_max = zeros(1, numberOfData);
gap_bits        = 0.5;
for i = 1:numberOfData
    mi = mutual_informations(i,:);
    index_time_almost_max = find(mi < mi(end) - gap_bits); % assume reaching equilibrium
    if isempty(index_time_almost_max)
        time_almost_max(i) = 1;
    else
        time_almost_max(i) = index_time_almost_max(end) + 1;
    end
end

figure()
mi = mutual_informations(:,end);
mi = mi(1:end-1);
plot(mi, -log10(time_almost_max(1:end-1)), '-o')
xlabel('MI')
ylabel('-log10(Time Required to Reach Final MI in days)')
title(sprintf('r0 = %0.2f, alpha = %0.2f, sigma = %0.2f, ext noise lvl = %0.2f', r0, alpha, sigmas(end), ext_noise_lvls(end)));
hold all;
plot(mutual_informations(end,end), -log10(time_almost_max(end)), '-ro', 'MarkerSize', 8, 'LineWidth', 8)
hold off;

figure()
plot(r0./Lmaxs(1:end-1), mi,'-o')
xlabel('R/L')
ylabel('MI')
title(sprintf('r0 = %0.2f, alpha = %0.2f, sigma = %0.2f, ext noise lvl = %0.2f', r0, alpha, sigmas(end), ext_noise_lvls(end)));
hold all;
plot(1/1000, mutual_informations(end,end), '-ro', 'MarkerSize', 8, 'LineWidth', 8)
hold off;

% % for i = 1:numberOfData
% %     diff_mi = diff(mutual_informations(i, :));
% %     sorted_diff_mi = sort(diff_mi, 'descend');
% %     learning_rates(i) = mean(sorted_diff_mi(1:2));
% % end
% figure()
% semilogx(r0./Lmaxs, learning_rates, '-o')
% xlabel('R/L')
% ylabel('Learning Rate')
% title(sprintf('r0 = %0.2f, alpha = %0.2f', r0, alpha));
% 
% figure()
% plot(mutual_informations(:, end), learning_rates, '-o')
% xlabel('MI')
% ylabel('Learning Rate')
% title(sprintf('r0 = %0.2f, alpha = %0.2f', r0, alpha));


% figure()
% for i = 1:length(sigmas)
%     plot(0:dt:finalDay, entropies(i, :), 'DisplayName', sprintf('sigma = %0.2f, ext noise lvls = %0.2f', sigmas(i), ext_noise_lvls(i)))
%     xlabel('time (days)')
%     ylabel('Entropy')
%     hold all;
% end
% hold off;
% title(sprintf('r0 = %0.2f, L = %0.2f, alpha = %0.2f', r0, Lmax, alpha));
% legend('-DynamicLegend')

