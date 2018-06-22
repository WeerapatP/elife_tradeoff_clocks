numberOfData = 12;
RoverLs      = 5*logspace(-1,1,numberOfData);

alpha = 5;
r0    = 0.5;
ensemblesize = 10000;
dt    = 0.001;
noise_correlation_time = 0.1;

% for external noise Appendix 5 Fig 1 d
sigma = 0;
ext_noise_lvl = 1;

% for internal noise Appendix 5 Fig 1 e
% sigma = 0.1;
% ext_noise_lvl = 0;

finalDay = 1000;

variance_gains  = zeros(1, numberOfData);
entropy_drops  = zeros(1, numberOfData);
sigma_es       = zeros(1, numberOfData);
mis            = zeros(1, numberOfData);
parfor i = 1:numberOfData
    i
    Lmax = r0/RoverLs(i);
    [mid_entropy, max_entropy, mis(i)] = simulate_attractor(alpha, r0, Lmax, ensemblesize, dt, noise_correlation_time, sigma, ext_noise_lvl, finalDay)
    variance_gains(i) = 2^(max_entropy) - 2^(mid_entropy);
    entropy_drops(i) = max_entropy - mid_entropy;
end

% Appendix 5 Fig 1
figure(102)
plot(log10(RoverLs), entropy_drops, '-o')
xlabel('R/L')
ylabel('Entropy Drop at Dusk')
set(gca, 'fontsize', 14)

% Fig 4 d
figure(103)
plot(log10(RoverLs), log10(variance_gains), '-o')
xlabel('R/L')
ylabel('log2(Variance) Gain during the day')
set(gca, 'fontsize', 14)

figure(104)
plot(RoverLs, mis, '-o')
xlabel('R/L')
ylabel('MI')
set(gca, 'fontsize', 14)

