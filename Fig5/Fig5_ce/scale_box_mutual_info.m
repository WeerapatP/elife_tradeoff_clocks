%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Matlab code to evolve phase
%  point(s) according to L-dependent
%  limit cycles. Internal noise are 
%  assumed to be gaussian. Langevine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Zhiyue Lu,
%  created on Jan 10, 2017
%  last modified on Jan. 20, 2017
%  zhiyuelu@gmail.com if any question
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delete(gcp('nocreate'))
load carbig

pi = 3.1415926535897932384626433;
pre_genesis=8000; %number of "days" of darkness for initialization
n_days=pre_genesis+8000;% Total number of days, including before genesis
hd=12;%12 hour daylight
frac_day=hd/24;
chain_length=24*n_days; % I assume everyday has 24 hours.
ensemblesize=40000; % number of independent uniformly generated points to run
%Lmax=1.0;   % this is the maximum intensity of light(maxshift of limit cycle)
maxLmax=Lmax;
Lmax1=3.0*Lmax;
L=0.0;
r0=1.5;   % this is the radius of the night limit cycle (non-shifted)
gamma=0.0; % this is the rate of ???
alpha=4.0;  % this is the rate of ???
%ctr={-r0-0.5*Lmax:0.02:r0+1.5*maxLmax -r0-1*maxLmax:0.02:r0+1*maxLmax};
ctrcoarse={-2*r0:0.05:2*r0+Lmax1 -Lmax1/2-2*r0:0.05:2*r0+Lmax1/2};
omega0=2*pi;  % this is omega of night, changing it doesn't change the fact
% that the period of a night limit cycle is 24 hours.
dt=0.0001;
sqrtdt=sqrt(dt);
sigma=sqrt(0.1);
xmax=0;
ymax=0;
xmin=0;
ymin=0;


x=zeros(ensemblesize,2);
xx=zeros(ensemblesize,1);
yy=zeros(ensemblesize,1);
dxx=zeros(ensemblesize,1);
dyy=zeros(ensemblesize,1);
%%%%%%%%set initial ensemble on the ring of darkness...
for i=1:ensemblesize
    x(i,1)=normrnd(0.0,sigma);
    x(i,2)=normrnd(0.0,sigma);
    radius=sqrt(x(i,1)^2+x(i,2)^2);
    x(i,1)=r0*x(i,1)/radius;
    x(i,2)=r0*x(i,2)/radius;
end
xx=x(:,1);
yy=x(:,2);


oldS=-10;
%pregenesis, preparation
for time=0:dt:pre_genesis
    
    L=0.0;
    dxx=(alpha*(1-sign(alpha)*((xx-L).^2+yy.^2)/r0^2).*(xx-L)-yy.*(omega0+gamma*((xx-L).^2+yy.^2-r0^2))).*dt + normrnd(0.0,sigma,ensemblesize,1)*sqrtdt;
    dyy=(alpha*(1-sign(alpha)*((xx-L).^2+yy.^2)/r0^2).*yy+(xx-L).*(omega0+gamma*((xx-L).^2+yy.^2-r0^2))).*dt + normrnd(0.0,sigma,ensemblesize,1)*sqrtdt;
    xx=xx+dxx;
    yy=yy+dyy;
    %the following part checks if the entropy at begning of this day is almost same with that at start of yesterday 
    if(mod(time,1)==0.0 & time>=5.0)
        N=hist3([xx,yy],ctrcoarse);
        N=N/ensemblesize;
        newS=-sum(N(N~=0).*log(N(N~=0)));
        if((abs(newS-oldS)/newS)<1E-3) break;
        else oldS=newS;
        end
    end
end
x=[xx,yy];




%initial state generated, below starts the day-night cycle
oldS=-100.0;
startrecord=0;
starttime=n_days+1;
emptyscal=0;






        

time



for time=pre_genesis:dt:n_days
    if (rem(time,1)<frac_day) L=Lmax;
    else L=0.0;
    end
    
    dxx=(alpha*(1-sign(alpha)*((xx-L).^2+yy.^2)/r0^2).*(xx-L)-yy.*(omega0+gamma*((xx-L).^2+yy.^2-r0^2))).*dt + normrnd(0.0,sigma,ensemblesize,1)*sqrtdt;
    dyy=(alpha*(1-sign(alpha)*((xx-L).^2+yy.^2)/r0^2).*yy+(xx-L).*(omega0+gamma*((xx-L).^2+yy.^2-r0^2))).*dt + normrnd(0.0,sigma,ensemblesize,1)*sqrtdt;
    xx=xx+dxx;
    yy=yy+dyy;

    %the following part checks if the entropy at begning of this day is almost same with that at start of yesterday 
    if(mod(time,1)==0.0 & (time-pre_genesis)>=5.0 & startrecord==0)
        N=hist3([xx,yy],ctrcoarse);
        N=N/ensemblesize;
        newS=-sum(N(N~=0).*log(N(N~=0)));
        if((abs(newS-oldS)/newS)<1E-3)  % this is when the condition is satisfied
            startrecord=1;
            starttime=time;
        else oldS=newS;
        end
    end
    
    %here we calculate 
    if(startrecord==1)
        xmax=max(xmax,max(xx));
        xmin=min(xmin,min(xx));
        ymax=max(ymax,max(yy));
        ymin=min(ymin,min(yy));
        % here we find the box boundary at a periodic state.
    end
    
    if(startrecord==1 & time-starttime>=(1.0-2*dt))
        break;
    end
end

ddx=0.05;
ddy=0.05;
ctrcoarse={xmin:ddx:xmax ymin:ddy:ymax};
N=hist3([xx,yy],ctrcoarse);

emptyN=zeros(size(N));

%for prints the shannon entropy_vs time
strr = sprintf('check_r0_%dalpha_%dsigma_%dLmax_%d.dat',r0,alpha,sigma,Lmax);
fileIDD = fopen(strr,'w');

%here is the last day for statistics:
for time=0:dt:2.0-2*dt
    if (rem(time,1)<frac_day) L=Lmax;
    else L=0.0;
    end
    
    dxx=(alpha*(1-sign(alpha)*((xx-L).^2+yy.^2)/r0^2).*(xx-L)-yy.*(omega0+gamma*((xx-L).^2+yy.^2-r0^2))).*dt + normrnd(0.0,sigma,ensemblesize,1)*sqrtdt;
    dyy=(alpha*(1-sign(alpha)*((xx-L).^2+yy.^2)/r0^2).*yy+(xx-L).*(omega0+gamma*((xx-L).^2+yy.^2-r0^2))).*dt + normrnd(0.0,sigma,ensemblesize,1)*sqrtdt;
    xx=xx+dxx;
    yy=yy+dyy;
    
    N=hist3([xx,yy],ctrcoarse);
    N=N/ensemblesize;
    shannon=-sum(N(N~=0).*log(N(N~=0)));
    largevar = max(eig(cov([xx yy])));
    
    fprintf(fileIDD,'%d\t%d\n',time,log(largevar));
    emptyscal=emptyscal+(shannon)*dt;
    emptyN=emptyN+N*dt/2.0;%%%%%%%%%%%!!!!!%%% careful here customized for two day final measurement
end

str=sprintf('two_dim_result.txt');
fileID = fopen(str,'a');
if (starttime>n_days) fprintf(fileID,'error\n');
end
sum(sum(emptyN));
mutualinformation=-emptyscal/2.0-sum(emptyN(emptyN~=0).*log(emptyN(emptyN~=0)));%%%%%%%%!!!!!%%%%%% careful here customized for two day final measurement
%fprintf(fileID,'# ro \t L \t alpha \t sigma \t mutualinformation\n');
fprintf(fileID,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',r0,Lmax,alpha,sigma,mutualinformation,emptyscal/2.0,xmax,xmin,ymax,ymin);
%%%%%%%%!!!!!%%%%%% careful here customized for two day final measurement




close all
clf

