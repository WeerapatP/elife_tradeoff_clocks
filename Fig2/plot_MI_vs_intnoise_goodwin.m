%% limit cycle
sim_length = 12;
Omegas   = logspace(2, 7, sim_length);
has_extnoise = true;
has_intnoise = true;

n = 9;
s2 =  108;
s1 = 92;
mis1 = zeros(1, sim_length);
signal_maxs  = s2*ones(1, sim_length);
signal_mins  = s1*ones(1, sim_length);
parfor i = 1:sim_length
    s_max = signal_maxs(i);
    s_min = signal_mins(i);
    Omega  = Omegas(i);
    mis1(i) = simulate_goodwin(s_min, s_max, n, Omega, has_extnoise, has_intnoise);
    ['sim i = ', num2str(i), ' is finished.']
end
mis1
legend('-DynamicLegend');
set(gca,'fontsize',18);

s2 =  120;
s1 = 80;
mis2 = zeros(1, sim_length);
signal_maxs  = s2*ones(1, sim_length);
signal_mins  = s1*ones(1, sim_length);
parfor i = 1:sim_length
    s_max = signal_maxs(i);
    s_min = signal_mins(i);
    Omega  = Omegas(i);
    mis2(i) = simulate_goodwin(s_min, s_max, n, Omega, has_extnoise, has_intnoise);
    ['sim i = ', num2str(i), ' is finished.']
end
mis2
legend('-DynamicLegend');
set(gca,'fontsize',18);

%% point attractor
n = 9;
s2 =  2.5;
s1 =  1;
mis3 = zeros(1, sim_length);
signal_maxs  = s2*ones(1, sim_length);
signal_mins  = s1*ones(1, sim_length);
parfor i = 1:sim_length
    s_max = signal_maxs(i);
    s_min = signal_mins(i);
    Omega  = Omegas(i);
    mis3(i) = simulate_goodwin(s_min, s_max, n, Omega, has_extnoise, has_intnoise);
    ['sim i = ', num2str(i), ' is finished.']
end
mis3
legend('-DynamicLegend');
set(gca,'fontsize',18);

%%
figure()
% plot(log10(Omegas), mis1, '-o', log10(Omegas), mis2, '-o', log10(Omegas), mis3, '-o')
plot(log10(Omegas), mis1, '-o', log10(Omegas), mis2, '-o',log10(Omegas), mis3, '-o')
xlabel('log10(Omegas)')
ylabel('MI')
legend('Cycle with small L', 'Cycle with large L', 'Point Attractor')
title('Goodwin')
set(gca,'fontsize',18);