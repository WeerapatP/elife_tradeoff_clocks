finalDay = 100;
dt    = 0.001;
gamma = 0;
ensemblesize = 100000;

% Limit Cycle & Point Attractor
alphas  = [  5      5    -5];
r0s     = [1        1  1000];
Lmaxs   = [0.1    0.5   1.5];

% Sorted of Fixed Variables
noise_correlation_time = 0.1;
ext_noise_lvl = 1;
sigmas = logspace(-2,0,12);
mis = zeros(length(alphas), length(sigmas));
        
parfor i = 1:length(alphas)
    for j = 1:length(sigmas)
        r0   = r0s(i);
        Lmax = Lmaxs(i);
        alpha = alphas(i);
        sigma = sigmas(j);
        mis(i, j) = simulate_attractor(alpha, r0, Lmax, ensemblesize, dt, noise_correlation_time, sigma, ext_noise_lvl, finalDay);
    end
end

%%
figure()
semilogx(sigmas, mis(1,:), '-o', sigmas, mis(2,:), '-o', sigmas, mis(3,:), '-o')

xlabel('\epsilon_{int}')
ylabel('MI')
legend('Cycle with small L', 'Cycle with large L', 'Point Attractor')