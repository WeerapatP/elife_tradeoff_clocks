A1=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_1.000000e-01.dat');
A6=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_1.100000e+00.dat');
A7=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_1.300000e+00.dat');
A8=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_1.500000e+00.dat');
A9=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_1.700000e+00.dat');
A10=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_1.900000e+00.dat');
A11=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_2.100000e+00.dat');
A12=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_2.300000e+00.dat');
A13=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_2.500000e+00.dat');
A14=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_2.700000e+00.dat');
A15=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_2.900000e+00.dat');
A2=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_3.000000e-01.dat');
A16=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_3.100000e+00.dat');
A3=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_5.000000e-01.dat');
A4=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_7.000000e-01.dat');
A5=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_9.000000e-01.dat');

A23=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_3.000000e-02.dat');
A22=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_3.333333e-02.dat');
A21=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_3.750000e-02.dat');
A20=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_4.285714e-02.dat');
A19=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_5.000000e-02.dat');
A18=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_6.000000e-02.dat');
A17=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_7.500000e-02.dat');
A25=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_1.500000e-01.dat');

A28=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_7.142857e-02.dat');
A27=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_8.823529e-02.dat');
A26=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_1.250000e-01.dat');
A24=importdata('check_r0_1.500000e+00alpha_4sigma_3.162278e-01Lmax_2.142857e-01.dat');


lab=size(A1);
lll=lab(1);
samplerate=199;

BLi=20:5:50;
BLmakeup = 1.5./BLi; 
Bx=[0.1:0.2:3.1 BLmakeup 1.5/7 1.5/10 1.5/12 1.5/17 1.5/21];
Bxx=1.5./Bx;

B(1)=max(A1(:,2))-min(A1(:,2));
B(2)=max(A2(:,2))-min(A2(:,2));
B(3)=max(A3(:,2))-min(A3(:,2));
B(4)=max(A4(:,2))-min(A4(:,2));
B(5)=max(A5(:,2))-min(A5(:,2));
B(6)=max(A6(:,2))-min(A6(:,2));
B(7)=max(A7(:,2))-min(A7(:,2));
B(8)=max(A8(:,2))-min(A8(:,2));
B(9)=max(A9(:,2))-min(A9(:,2));
B(10)=max(A10(:,2))-min(A10(:,2));
B(11)=max(A11(:,2))-min(A11(:,2));
B(12)=max(A12(:,2))-min(A12(:,2));
B(13)=max(A13(:,2))-min(A13(:,2));
B(14)=max(A14(:,2))-min(A14(:,2));
B(15)=max(A15(:,2))-min(A15(:,2));
B(16)=max(A16(:,2))-min(A16(:,2));
B(17)=max(A17(:,2))-min(A17(:,2));
B(18)=max(A18(:,2))-min(A18(:,2));
B(19)=max(A19(:,2))-min(A19(:,2));
B(20)=max(A20(:,2))-min(A20(:,2));
B(21)=max(A21(:,2))-min(A21(:,2));
B(22)=max(A22(:,2))-min(A22(:,2));
B(23)=max(A23(:,2))-min(A23(:,2));
B(24)=max(A24(:,2))-min(A24(:,2));
B(25)=max(A25(:,2))-min(A25(:,2));
B(26)=max(A26(:,2))-min(A26(:,2));
B(27)=max(A27(:,2))-min(A27(:,2));
B(28)=max(A28(:,2))-min(A28(:,2));

figure(1)
plot([Bxx(28) Bxx(17) Bxx(27) Bxx(1) Bxx(26:-1:24) Bxx(2:16)],[B(28) B(17) B(27) B(1) B(26:-1:24) B(2:16)]./log(2),'x-')
xlabel('R/L')
ylabel('Shannon Entropy Drop (bit)')




clear B;

B(1)=mean(A1(:,2));
B(2)=mean(A2(:,2));
B(3)=mean(A3(:,2));
B(4)=mean(A4(:,2));
B(5)=mean(A5(:,2));
B(6)=mean(A6(:,2));
B(7)=mean(A7(:,2));
B(8)=mean(A8(:,2));
B(9)=mean(A9(:,2));
B(10)=mean(A10(:,2));
B(11)=mean(A11(:,2));
B(12)=mean(A12(:,2));
B(13)=mean(A13(:,2));
B(14)=mean(A14(:,2));
B(15)=mean(A15(:,2));
B(16)=mean(A16(:,2));
B(17)=mean(A17(:,2));
B(18)=mean(A18(:,2));
B(19)=mean(A19(:,2));
B(20)=mean(A20(:,2));
B(21)=mean(A21(:,2));
B(22)=mean(A22(:,2));
B(23)=mean(A23(:,2));
B(24)=mean(A24(:,2));
B(25)=mean(A25(:,2));
B(26)=mean(A26(:,2));
B(27)=mean(A27(:,2));
B(28)=mean(A28(:,2));

B(1)=mean(A1(:,2));
B(2)=mean(A2(:,2));
B(3)=mean(A3(:,2));
B(4)=mean(A4(:,2));
B(5)=mean(A5(:,2));
B(6)=mean(A6(:,2));
B(7)=mean(A7(:,2));
B(8)=mean(A8(:,2));
B(9)=mean(A9(:,2));
B(10)=mean(A10(:,2));
B(11)=mean(A11(:,2));
B(12)=mean(A12(:,2));
B(13)=mean(A13(:,2));
B(14)=mean(A14(:,2));
B(15)=mean(A15(:,2));
B(16)=mean(A16(:,2));
B(17)=mean(A17(:,2));
B(18)=mean(A18(:,2));
B(19)=mean(A19(:,2));
B(20)=mean(A20(:,2));
B(21)=mean(A21(:,2));
B(22)=mean(A22(:,2));
B(23)=mean(A23(:,2));
B(24)=mean(A24(:,2));
B(25)=mean(A25(:,2));
B(26)=mean(A26(:,2));
B(27)=mean(A27(:,2));
B(28)=mean(A28(:,2));


figure(2)
plot([Bxx(28) Bxx(17) Bxx(27) Bxx(1) Bxx(26:-1:24) Bxx(2:16)],[B(28) B(17) B(27) B(1) B(26:-1:24) B(2:16)].*log(2),'x-')
xlabel('R/L')
ylabel('mean Shannon entropy (bit)')






