% simulate_millar(0.02,0.02,0.26,100,false,false)
% simulate_millar(0.005,0.005,0.26,100,false,false)
% simulate_millar(0,0,0.26,100,false,false)
% simulate_millar(-0.01,-0.01,0.26,100,false,false)

has_extnoise = true;
has_intnoise = true;
sim_length = 12;
Omegas   = logspace(2, 7, sim_length);

signal_min = 0;
% signal_max = 0.007;
signal_max = 0.01;
mis1 = zeros(1, sim_length);
v1 = 0.26;
parfor i = 1:sim_length
    Omega = Omegas(i);
    mis1(i)  = simulate_millar(signal_min, signal_max, v1, Omega, has_extnoise, has_intnoise);
end
mis1

% signal_min = -0.005;
signal_min = 0;
% signal_max = 0.02;
signal_max = 0.05;
mis2 = zeros(1, sim_length);
v1 = 0.26;
parfor i = 1:sim_length
    Omega = Omegas(i);
    mis2(i)  = simulate_millar(signal_min, signal_max, v1, Omega, has_extnoise, has_intnoise);
end
mis2

signal_min = 0;
signal_max = 0.2;
mis3 = zeros(1, sim_length);
v1 = 0.05;
parfor i = 1:sim_length
    Omega = Omegas(i);
    mis3(i)  = simulate_millar(signal_min, signal_max, v1, Omega, has_extnoise, has_intnoise);
end
mis3

figure()
plot(log10(Omegas), mis1, '-o', log10(Omegas), mis2, '-o', log10(Omegas), mis3, '-o')
xlabel('log10(Omegas)')
ylabel('MI')
legend('Cycle with small L', 'Cycle with large L', 'Point Attractor')
set(gca,'fontsize',18);
title('Millar');