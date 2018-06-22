%%
has_extnoise = true;
has_intnoise = true;
sim_length = 12;
Omegas   = logspace(2, 10, sim_length);

s2 = 0.0106;
s1 = 0.0105;
mis1 = zeros(1, sim_length);
parfor i = 1:sim_length
    Omega  = Omegas(i);
    mis1(i) = simulate_cellcycle(s1, s2, Omega, has_extnoise, has_intnoise);
    ['sim i = ', num2str(i), ' is finished.']
end
mis1

s2 = 0.0111;
s1 = 0.0105;
mis2 = zeros(1, sim_length);
parfor i = 1:sim_length
    Omega  = Omegas(i);
    mis2(i) = simulate_cellcycle(s1, s2, Omega, has_extnoise, has_intnoise);
    ['sim i = ', num2str(i), ' is finished.']
end
mis2

s2 = 0.009;
s1 = 0.00;
mis3 = zeros(1, sim_length);
parfor i = 1:sim_length
    Omega  = Omegas(i);
    mis3(i) = simulate_cellcycle(s1, s2, Omega, has_extnoise, has_intnoise);
    ['sim i = ', num2str(i), ' is finished.']
end
mis3

figure()
plot(log10(Omegas), mis1, '-o', log10(Omegas), mis2, '-o', log10(Omegas), mis3, '-o')
xlabel('log10(Omegas)')
ylabel('MI')
title('Cell Cycle')
legend('Cycle Small', 'Cycle Big', 'Point Attractor')
set(gca,'fontsize',18);