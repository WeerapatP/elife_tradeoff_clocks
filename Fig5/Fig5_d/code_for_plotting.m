C7=importdata('sigma_1.600000e+00_alpha_-1.251181e+00.dat');
C6=importdata('sigma_1.600000e+00_alpha_-1.565454e+00.dat');
C5=importdata('sigma_1.600000e+00_alpha_-1.958667e+00.dat');
C8=importdata('sigma_1.600000e+00_alpha_-1.dat');
C4=importdata('sigma_1.600000e+00_alpha_-2.450647e+00.dat');
C3=importdata('sigma_1.600000e+00_alpha_-3.066203e+00.dat');
C2=importdata('sigma_1.600000e+00_alpha_-3.836375e+00.dat');
C1=importdata('sigma_1.600000e+00_alpha_-4.800000e+00.dat');
A7=importdata('sigma_4.000000e-01_alpha_-1.251181e+00.dat');
A6=importdata('sigma_4.000000e-01_alpha_-1.565454e+00.dat');
A5=importdata('sigma_4.000000e-01_alpha_-1.958667e+00.dat');
A8=importdata('sigma_4.000000e-01_alpha_-1.dat');
A4=importdata('sigma_4.000000e-01_alpha_-2.450647e+00.dat');
A3=importdata('sigma_4.000000e-01_alpha_-3.066203e+00.dat');
A2=importdata('sigma_4.000000e-01_alpha_-3.836375e+00.dat');
A1=importdata('sigma_4.000000e-01_alpha_-4.800000e+00.dat');
B7=importdata('sigma_8.000000e-01_alpha_-1.251181e+00.dat');
B6=importdata('sigma_8.000000e-01_alpha_-1.565454e+00.dat');
B5=importdata('sigma_8.000000e-01_alpha_-1.958667e+00.dat');
B8=importdata('sigma_8.000000e-01_alpha_-1.dat');
B4=importdata('sigma_8.000000e-01_alpha_-2.450647e+00.dat');
B3=importdata('sigma_8.000000e-01_alpha_-3.066203e+00.dat');
B2=importdata('sigma_8.000000e-01_alpha_-3.836375e+00.dat');
B1=importdata('sigma_8.000000e-01_alpha_-4.800000e+00.dat');



MA(1)=mean(A1(:,2));
MA(2)=mean(A2(:,2));
MA(3)=mean(A3(:,2));
MA(4)=mean(A4(:,2));
MA(5)=mean(A5(:,2));
MA(6)=mean(A6(:,2));
MA(7)=mean(A7(:,2));
MA(8)=mean(A8(:,2));


MB(1)=mean(B1(:,2));
MB(2)=mean(B2(:,2));
MB(3)=mean(B3(:,2));
MB(4)=mean(B4(:,2));
MB(5)=mean(B5(:,2));
MB(6)=mean(B6(:,2));
MB(7)=mean(B7(:,2));
MB(8)=mean(B8(:,2));


MC(1)=mean(C1(:,2));
MC(2)=mean(C2(:,2));
MC(3)=mean(C3(:,2));
MC(4)=mean(C4(:,2));
MC(5)=mean(C5(:,2));
MC(6)=mean(C6(:,2));
MC(7)=mean(C7(:,2));
MC(8)=mean(C8(:,2));




Alpha=-logspace(-log10(5/24),-log10(24/24),8);
tauxx = -24./Alpha;

%Bxx=-0.4*2^jj;
figure(5)
semilogx(tauxx,log2(MA)-2*log2(0.4),'o-',tauxx,log2(MB)-2*log2(0.8),'*-',tauxx,log2(MC)-2*log2(1.6),'x-')
xlabel('tau/hr')
ylabel('log of mean variance -log2(D) (bit) ')
legend('D=0.16','D=0.64', 'D=2.56')










