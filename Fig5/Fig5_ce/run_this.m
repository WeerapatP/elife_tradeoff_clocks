clear all
clc
clf
close all



for jj=20:5:50
    Lmax=1.5/jj;
    run('scale_box_mutual_info.m')
    clear emptyN;
    clear N;
end

   

clear all
clc
clf
close all
 for jj=1:1:16
     Lmax=0.2*jj-0.1;
     run('scale_box_mutual_info.m')
     clear emptyN;
     clear N;
 end
% 
%     Lmax=6;
%     run('scale_box_mutual_info.m')
%     clear emptyN;
%     clear N;


    Lmax=1.5/7;
    run('scale_box_mutual_info.m')
    clear emptyN;
    clear N;
    
    Lmax=1.5/12;
    run('scale_box_mutual_info.m')
    clear emptyN;
    clear N;
    
    Lmax=1.5/17;
    run('scale_box_mutual_info.m')
    clear emptyN;
    clear N;
    
    Lmax=1.5/21;
    run('scale_box_mutual_info.m')
    clear emptyN;
    clear N;