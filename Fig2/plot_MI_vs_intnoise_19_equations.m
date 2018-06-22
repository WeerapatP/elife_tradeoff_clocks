% simulate_19_equations(1.1,1.1,100,false,false)
% simulate_19_equations(1,1,100,false,false)
% simulate_19_equations(0.3,1.0,100,false,false)
% simulate_19_equations(0.1,1.0,100,false,false)
% simulate_19_equations(0.1,0.9,100,false,false)
% simulate_19_equations(1.07, 1.09, 100, false, false, false)
% simulate_19_equations(1.07, 1.15, 100, false, false, false)
% simulate_19_equations(0, 1.5, 100, false, false, true)
% 

%%
has_extnoise = true;
has_intnoise = true;
sim_length = 12;
Omegas   = logspace(1, 7, sim_length);

s2 = 1.09;
s1 = 1.07;
n=3;
mis1 = zeros(1, sim_length);
parfor i = 1:sim_length
    Omega  = Omegas(i);
    mis1(i) = simulate_19_equations(s1, s2, Omega, has_extnoise, has_intnoise, false);
    ['sim i = ', num2str(i), ' is finished.']
end
mis1

s2 = 1.15;
s1 = 1.07;
mis2 = zeros(1, sim_length);
parfor i = 1:sim_length
    Omega  = Omegas(i);
    mis2(i) = simulate_19_equations(s1, s2, Omega, has_extnoise, has_intnoise, false);
    ['sim i = ', num2str(i), ' is finished.']
end
mis2

s2 = 1.5;
s1 = 0;
mis3 = zeros(1, sim_length);
parfor i = 1:sim_length
    Omega  = Omegas(i);
    mis3(i) = simulate_19_equations(s1, s2, Omega, has_extnoise, has_intnoise, true);
    ['sim i = ', num2str(i), ' is finished.']
end
mis3

figure()
plot(log10(Omegas), mis1, '-o', log10(Omegas), mis2, '-o', log10(Omegas), mis3, '-o')
xlabel('log10(Omegas)')
ylabel('MI')
title('Mammalian Clock')
legend('Cycle Small', 'Cycle Big', 'Point Attractor')
set(gca,'fontsize',18);