clear
full_pdf_name = './c.pdf';
figure_dest = './PDF/';

%openfig('Corrected Fig 2 d entropy drop_vs_RoverL.fig');
openfig('priliminary_c.fig');
h2 = findobj(gca,'Type','line');
roverl = get(h2,'Xdata');
deltaEntropy = get(h2,'Ydata');
%
close all;


% Plot them together
fig=figure(2); 
clf;

% Set size of fonts
xtickSize = 16;  ytickSize = 16;
xlabelSize = 20; ylabelSize = 20;

%xticks([0 0.25 0.5 0.75 1.0])
%xticklabels({'6','12','18','24','6'})
%set(gca,'Position',[0.1 0.55 0.65 0.35])


% Mathematica data
s2minus1TheoryFromMathematica = [[60., 37.5, 30., 12., 6., 3, 2., 3/2, 1.2, 1, 0.857143, 0.8, 3/4, ...
  0.666667, 3/5, 1/2, 3/8]; [0.0168056, 0.0270222, 0.0338889, ...
  0.0868056, 0.180556, 0.388889, 0.625, 0.888889, 1.18056, 1.5, ...
  1.84722, 2.03125, 2.22222, 2.625, 3.05555, 4., 6.22222]];

brownrgb=[139,69,19]/256;

loglog(s2minus1TheoryFromMathematica(1,:),s2minus1TheoryFromMathematica(2,:),'-','Color',brownrgb,'LineWidth',0.5)
hold all;

loglog(roverl,1.2*((2.^(deltaEntropy))-1),'xg','LineWidth',0.5,'MarkerSize',12)
hold off

xlabel('R/L');
ylabel('(s^2-1) and  e*(2\^ (zhiuye data)-1)');
%
pbaspect([1 1.2 1]);

xlim([0.4 32])
ylim([0.02 3])
%
xlabel('R/L');
set(get(gca,'XLabel'), 'FontSize', xlabelSize);
set(gca,'XTick',[0.5 1 2 4 8 16 32 ]);
%set(gca,'XTickLabel',{'0.5','1','5','15'},'FontSize', xtickSize);
set(gca,'YTick',[0.01 0.02 0.04 0.08 0.16 0.32 0.64 1.28 2.56]);
ylabel('Dawn variance drop s^2 - 1');
set(get(gca,'YAxis'),'FontSize', ytickSize)
set(get(gca,'YLabel'), 'FontSize', ylabelSize);
set(get(gca,'XAxis'),'FontSize', xtickSize)
set(get(gca,'XLabel'), 'FontSize', xlabelSize);

 saveas(fig, full_pdf_name);
% close all;


%% New panel f:
%
%
%
%

clear all
full_pdf_name = './e.pdf';
figure_dest = './PDF/';

openfig('priliminary_e.fig');
h4 = findobj(gca,'Type','line');
x = get(h4,'Xdata');
y = get(h4,'Ydata');
close all;


% Plot them together
fig=figure(4); 
clf;

darkgreen=[0.1 1 0.1]/2.2;
%plot(x{1,1},y{1,1},'*r','LineWidth',0.5,'MarkerSize',14);
%hold on
%plot(x,-y,'x','Color',darkgreen,'LineWidth',0.5);
%hold off
%
% Set size of fonts
xtickSize = 14;  ytickSize = 14;
xlabelSize = 16; ylabelSize = 16;

%xticks([0 0.25 0.5 0.75 1.0]
%xticklabels({'6','12','18','24','6'})
%set(gca,'Position',[0.1 0.55 0.65 0.35])


%

% Mathematica data
s2minus1TheoryFromMathematica = [[60., 37.5, 30., 12., 6., 3, 2., 3/2, 1.2, 1, 0.857143, 0.8, 3/4, ...
  0.666667, 3/5, 1/2, 3/8]; [0.0168056, 0.0270222, 0.0338889, ...
  0.0868056, 0.180556, 0.388889, 0.625, 0.888889, 1.18056, 1.5, ...
  1.84722, 2.03125, 2.22222, 2.625, 3.05555, 4., 6.22222]];

roverl = s2minus1TheoryFromMathematica(1,:);
s2 = s2minus1TheoryFromMathematica(2,:)+1;
varss = (s2+1)./(s2.^2-1);
clf;
% Superpose theory
hold all;
plot((roverl(3:end-4)),y(3)/log(2)^2+5.2+1*log2(1./varss(3:end-4)),'-k');

plot([x 10],(y(3)+[-y 1.3292*log(2)])/log(2)^2,'x','Color',darkgreen,'LineWidth',0.5);


plot(0,2.344/log(2)^2,'*r','LineWidth',0.5,'MarkerSize',14);


xlim([-1 22]);
ylim([-0.5 5.5]);
pbaspect([1 1.1 1]);


xlabel('R/L');
%
set(get(gca,'XLabel'), 'FontSize', xlabelSize);
set(get(gca,'XAxis'),'FontSize', xtickSize)


%set(gca,'XTick',[0:5:15]);
%set(gca,'XTickLabel',{'0','5','10','15'},'FontSize', xtickSize);
%set(gca,'YTick',[0:.5:2.5]);
ylabel('Precision (bits)');
%lgd=legend('D=0.1','D=25.6');
lgd.FontSize = 15;
set(get(gca,'YAxis'),'FontSize', ytickSize)
set(get(gca,'YLabel'), 'FontSize', ylabelSize);

saveas(fig, full_pdf_name);

%% compare theory+data on log scale to debug:
% Plot theory
figure;
clf; hold all;
plot(log(roverl),5+1*log2(1./varss),'-k');
% plot data:
plot([log(x) log(10)],[-y/(log(2))^2 1.3292/(log(2))],'O','LineWidth',0.5);
xlabel('log (r/l)');
ylabel(' precision (bits)');
xlim([0 4])
%%