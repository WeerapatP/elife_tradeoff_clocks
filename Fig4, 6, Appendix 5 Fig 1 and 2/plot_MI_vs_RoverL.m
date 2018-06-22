alpha = 5;
r0    = 0.5;
ensemblesize = 10000;
dt    = 0.001;
noise_correlation_time = 0.1;
ext_noise_lvl = 1;
sigma = 0;
finalDay = 200;
numberOfData = 8;
RoverLs      = logspace(0, 1, numberOfData);
mis          = zeros(1, numberOfData);

parfor i = 1:numberOfData
    'data'
    i
    Lmax = r0/RoverLs(i);
    [~,~,~,mis(i)] = simulate_attractor(alpha, r0, Lmax, ensemblesize, dt, noise_correlation_time, sigma, ext_noise_lvl, finalDay)
end

figure()
semilogx(RoverLs, mis, '-o')
xlabel('R/L')
ylabel('MI')
title(['r0 = ', num2str(r0), ', alpha = ', num2str(alpha), ', sigma = ', num2str(sigma)])


