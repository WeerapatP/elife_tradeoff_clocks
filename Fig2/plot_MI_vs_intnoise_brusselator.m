has_extnoise = true;
has_intnoise = true;
sim_length = 12;
Omegas   = logspace(2, 9, sim_length);

signal_min = 2.2;
signal_max = 2.25;
mis1 = zeros(1, sim_length);
parfor i = 1:sim_length
    Omega = Omegas(i);
    mis1(i)  = simulate_brusselator(signal_min, signal_max, Omega, has_extnoise, has_intnoise);
end
mis1

signal_min = 2.2;
signal_max = 2.8;
mis2 = zeros(1, sim_length);
parfor i = 1:sim_length
    Omega = Omegas(i);
    mis2(i)  = simulate_brusselator(signal_min, signal_max, Omega, has_extnoise, has_intnoise);
end
mis2

signal_min = 0.5;
signal_max = 1.8;
mis3 = zeros(1, sim_length);
parfor i = 1:sim_length
    Omega = Omegas(i);
    mis3(i)  = simulate_brusselator(signal_min, signal_max, Omega, has_extnoise, has_intnoise);
end
mis3

figure()
plot(log10(Omegas), mis1, '-o', log10(Omegas), mis2, '-o', log10(Omegas), mis3, '-o')
xlabel('log10(Omegas)')
ylabel('MI')
legend('Cycle with small L', 'Cycle with large L', 'Point Attractor')
set(gca,'fontsize',18);
title('Brusselator');