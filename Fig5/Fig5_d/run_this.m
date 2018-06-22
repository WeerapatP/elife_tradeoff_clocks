clear all
clc
clf
close all
sigma=sqrt(0.16);
Alphaxx=-logspace(-log10(5/24),-log10(24/24),8);
for jj=1:1:8
    alpha=Alphaxx(jj);
    run('scale_box_mutual_info.m')
    clear emptyN;
    clear N;
end

clear all
clc
clf
close all
sigma=sqrt(0.64);
Alphaxx=-logspace(-log10(5/24),-log10(24/24),8);
for jj=1:1:8
    alpha=Alphaxx(jj);
    run('scale_box_mutual_info.m')
    clear emptyN;
    clear N;
end

clear all
clc
clf
close all
sigma=sqrt(2.56);
Alphaxx=-logspace(-log10(5/24),-log10(24/24),8);
for jj=1:1:8
    alpha=Alphaxx(jj);
    run('scale_box_mutual_info.m')
    clear emptyN;
    clear N;
end


