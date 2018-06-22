% 0.65 - 17.9400
% 0.7 - 18.06
% 0.69 - 18.0667
% simulate_goldbeter(0.69,0.71,1000,true,true) %0.4668
% simulate_goldbeter(0.65,0.75,1000,true,true) %1.9605
% simulate_goldbeter(0.2,0.55,1000,true,true) % 1.7782

% simulate_goldbeter(0.69,0.71,100,true,true) %0.0308
% simulate_goldbeter(0.65,0.75,100,true,true) %0.29
% simulate_goldbeter(0.2,0.55,100,true,true) %1.0108

%% limit cycle
has_extnoise = true;
has_intnoise = true;
sim_length = 12;
Omegas   = logspace(1, 7, sim_length);

v1 = 0.695;
v2 = 0.705;
mis1 = zeros(1, sim_length);
parfor(i = 1:sim_length)
    Omega  = Omegas(i);
    mis1(i) = simulate_goldbeter(v1, v2, Omega, has_extnoise, has_intnoise);
    ['sim i = ', num2str(i), ' is finished.']
end
mis1
% 0.695 0.705
% mis1 = 0.0320    0.0825    0.3471    1.1216    2.0084    2.7015    3.2409    3.5714    3.7591    3.8448    3.8891    3.9055


v1 = 0.6;
v2 = 0.9;
mis2 = zeros(1, sim_length);
parfor(i = 1:sim_length)
    Omega  = Omegas(i);
    mis2(i) = simulate_goldbeter(v1, v2, Omega, has_extnoise, has_intnoise);
    ['sim i = ', num2str(i), ' is finished.']
end
mis2
% 0.6 0.9
% mis2 = 1.0539    1.8728    2.4957    2.9681    3.2657    3.4158    3.4800    3.5097    3.5075    3.5230    3.5208    3.5198
% 0.6 0.8
% mis2 = 0.7459    1.6035    2.3725    2.9164    3.2869    3.4907    3.5841    3.6107    3.6284    3.6347    3.6377    3.6389

%% point attractor
v2 = 0.55;
% v1 = 0.2;
v1 = 0.05;
mis3 = zeros(1, sim_length);
parfor(i = 1:sim_length)
    Omega  = Omegas(i);
    mis3(i) = simulate_goldbeter(v1, v2, Omega, has_extnoise, has_intnoise);
    ['sim i = ', num2str(i), ' is finished.']
end
mis3
% 0.2 0.55 - 1.0007    1.4136    1.7282    1.9304    2.0539    2.1193    2.1372    2.1618    2.1613    2.1708    2.1628    2.1622
%%
figure()
plot(log10(Omegas), mis1, '-o', log10(Omegas), mis2, '-o', log10(Omegas), mis3, '-o')
xlabel('log10(Omegas)')
ylabel('MI')
legend('Cycle with small L - vs = (0.69, 0.71)', 'Cycle with large L - vs = (0.6, 0.9)', 'Point attractor - vs = (0.05, 0.55)')
title('Goldbeter')
set(gca,'fontsize',18);