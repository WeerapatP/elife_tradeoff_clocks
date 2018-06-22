has_extnoise = true;
has_intnoise = true;
sim_length = 12;
Omegas   = logspace(1, 6, sim_length);

s2 = 5.9/5.2;
s1 = 4.5/5.2;
n=3;
alpha = 5.2;
mis1 = zeros(1, sim_length);
parfor i = 1:sim_length
    Omega  = Omegas(i);
    mis1(i) = simulate_repressilator(s1, s2, alpha, n, Omega, has_extnoise, has_intnoise);
    ['sim i = ', num2str(i), ' is finished.']
end
mis1

s2 = 10/5.2;
s1 = 3.5/5.2;
n=3;
alpha = 5.2;
mis2 = zeros(1, sim_length);
parfor i = 1:sim_length
    Omega  = Omegas(i);
    mis2(i) = simulate_repressilator(s1, s2, alpha, n, Omega, has_extnoise, has_intnoise);
    ['sim i = ', num2str(i), ' is finished.']
end
mis2

s2 = 2.4/1.9;
s1 = 0;
n=3;
alpha = 1.9;
mis3 = zeros(1, sim_length);
parfor i = 1:sim_length
    Omega  = Omegas(i);
    mis3(i) = simulate_repressilator(s1, s2, alpha, n, Omega, has_extnoise, has_intnoise);
    ['sim i = ', num2str(i), ' is finished.']
end
mis3

figure()
plot(log10(Omegas), mis1, '-o', log10(Omegas), mis2, '-o', log10(Omegas), mis3, '-o')
xlabel('log10(Omegas)')
ylabel('MI')
title('Repressilator')
legend('Cycle Small', 'Cycle Big', 'Point Attractor')
set(gca,'fontsize',18);
% signal2 = 0.1;
% signal1 = -0.1;
% n=2;
% alpha = 10;
% sigmal_maxs  = signal2*ones(1, sim_length);
% signal_mins  = signal1*ones(1, sim_length);
% mis1 = zeros(1, sim_length);
% parfor i = 1:sim_length
%     signal_max = sigmal_maxs(i);
%     signal_min = signal_mins(i);
%     Omega  = Omegas(i);
%     mis1(i) = simulate_repressilator(signal_min, signal_max, n, Omega, has_extnoise, has_intnoise);
%     ['sim i = ', num2str(i), ' is finished.']
% end
% mis1
% 
% signal2 = 0.7;
% signal1 = -0.3;
% n=3;
% sigmal_maxs  = signal2*ones(1, sim_length);
% signal_mins  = signal1*ones(1, sim_length);
% mis2 = zeros(1, sim_length);
% parfor i = 1:sim_length
%     signal_max = sigmal_maxs(i);
%     signal_min = signal_mins(i);
%     Omega  = Omegas(i);
%     mis2(i) = simulate_repressilator(signal_min, signal_max, n, Omega, has_extnoise, has_intnoise);
%     ['sim i = ', num2str(i), ' is finished.']
% end
% mis2