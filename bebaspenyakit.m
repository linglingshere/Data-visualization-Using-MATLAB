% PROGRAM BEBAS PENYAKIT (R0<1)
clc
clear
format long

tend=30;
h=0.001;
t=0:h:tend;
N=length(t);

%Nilai Parameter
L = 10; %lambda
b = 0.2; %beta
ts = 0.7; %theta 1
td = 0.1; %theta 2
tt = 0.3; %theta 3
mu = 0.3; %mu
d = 0.3; %delta
g = 0.2; %%gamma
s = 0.2; %sigma

%Vektor nol
S=0*t;
E=0*t;
I=0*t;
R=0*t;
 
%Nilai awal
S(1) = 100; 
E(1) = 15; 
I(1) = 20; 
R(1) = 4;
 
%Mencari 2 nilai awal dengan Metode Runge-Kutta-Fehlberg (RKF45) orde 5 
for i=1: N-1
    k1S = h*(L-b*(1-ts)*S(i)*(I(i)+td*E(i))-mu*S(i));
    k1E = h*(b*(1-ts)*S(i)*(I(i)+td*E(i))-(d+mu)*E(i));
    k1I = h*(d*E(i)-(g+s+mu)*I(i));
    k1R = h*((g+tt*s)*I(i)-mu*R(i));

    k2S = h*(L-b*(1-ts)*(S(i)+k1S/4)*((I(i)+k1I/4)+td*(E(i)+k1E/4))-mu*(S(i)+k1S/4));
    k2E = h*(b*(1-ts)*(S(i)+k1S/4)*((I(i)+k1I/4)+td*(E(i)+k1E/4))-(d+mu)*(E(i)+k1E/4));
    k2I = h*(d*(E(i)+k1E/4)-(g+s+mu)*(I(i)+k1I/4));
    k2R = h*((g+tt*s)*(I(i)+k1I/4)-mu*(R(i)+k1R/4));

    k3S = h*(L-b*(1-ts)*(S(i)+(3/32)*k1S+(9/32)*k2S)*((I(i)+(3/32)*k1I+(9/32)*k2I)+td*(E(i)+(3/32)*k1E+(9/32)*k2E))-mu*(S(i)+(3/32)*k1S+(9/32)*k2S));
    k3E = h*(b*(1-ts)*(S(i)+(3/32)*k1S+(9/32)*k2S)*((I(i)+(3/32)*k1I+(9/32)*k2I)+td*(E(i)+(3/32)*k1E+(9/32)*k2E))-(d+mu)*(E(i)+(3/32)*k1E+(9/32)*k2E));
    k3I = h*(d*(E(i)+(3/32)*k1E+(9/32)*k2E)-(g+s+mu)*(I(i)+(3/32)*k1I+(9/32)*k2I));
    k3R = h*((g+tt*s)*(I(i)+(3/32)*k1I+(9/32)*k2I)-mu*(R(i)+(3/32)*k1R+(9/32)*k2R));

    k4S = h*(L-b*(1-ts)*(S(i)+(1932/2197)*k1S-(7200/2197)*k2S+(7296/2197)*k3S)*((I(i)+(1932/2197)*k1I-(7200/2197)*k2I+(7296/2197)*k3I)+td*(E(i)+(1932/2197)*k1E-(7200/2197)*k2E+(7296/2197)*k3E))-mu*(S(i)+(1932/2197)*k1S-(7200/2197)*k2S+(7296/2197)*k3S));
    k4E = h*(b*(1-ts)*(S(i)+(1932/2197)*k1S-(7200/2197)*k2S+(7296/2197)*k3S)*((I(i)+(1932/2197)*k1I-(7200/2197)*k2I+(7296/2197)*k3I)+td*(E(i)+(1932/2197)*k1E-(7200/2197)*k2E+(7296/2197)*k3E))-(d+mu)*(E(i)+(1932/2197)*k1E-(7200/2197)*k2E+(7296/2197)*k3E));
    k4I = h*(d*(E(i)+(1932/2197)*k1E-(7200/2197)*k2E+(7296/2197)*k3E)-(g+s+mu)*(I(i)+(1932/2197)*k1I-(7200/2197)*k2I+(7296/2197)*k3I));
    k4R = h*((g+tt*s)*(I(i)+(1932/2197)*k1I-(7200/2197)*k2I+(7296/2197)*k3I)-mu*(R(i)+(1932/2197)*k1R-(7200/2197)*k2R+(7296/2197)*k3R));

    k5S = h*(L-b*(1-ts)*(S(i)+(439/216)*k1S-8*k2S+(3680/513)*k3S-(845/4104)*k4S)*((I(i)+(439/216)*k1I-8*k2I+(3680/513)*k3I-(845/4104)*k4I)+td*(E(i)+(439/216)*k1E-8*k2E+(3680/513)*k3E-(845/4104)*k4E))-mu*(S(i)+(439/216)*k1S-8*k2S+(3680/513)*k3S-(845/4104)*k4S));
    k5E = h*(b*(1-ts)*(S(i)+(439/216)*k1S-8*k2S+(3680/513)*k3S-(845/4104)*k4S)*((I(i)+(439/216)*k1I-8*k2I+(3680/513)*k3I-(845/4104)*k4I)+td*(E(i)+(439/216)*k1E-8*k2E+(3680/513)*k3E-(845/4104)*k4E))-(d+mu)*(E(i)+(439/216)*k1E-8*k2E+(3680/513)*k3E-(845/4104)*k4E));
    k5I = h*(d*(E(i)+(439/216)*k1E-8*k2E+(3680/513)*k3E-(845/4104)*k4E)-(g+s+mu)*(I(i)+(439/216)*k1I-8*k2I+(3680/513)*k3I-(845/4104)*k4I));
    k5R = h*((g+tt*s)*(I(i)+(439/216)*k1I-8*k2I+(3680/513)*k3I-(845/4104)*k4I)-mu*(R(i)+(439/216)*k1R-8*k2R+(3680/513)*k3R-(845/4104)*k4R));

    k6S = h*(L-b*(1-ts)*(S(i)-(8/27)*k1S+2*k2S-(3544/2565)*k3S+(1859/4104)*k4S-(11/40)*k5S)*((I(i)-(8/27)*k1I+2*k2I-(3544/2565)*k3I+(1859/4104)*k4I-(11/40)*k5I)+td*(E(i)-(8/27)*k1E+2*k2E-(3544/2565)*k3E+(1859/4104)*k4E-(11/40)*k5E))-mu*(S(i)-(8/27)*k1S+2*k2S-(3544/2565)*k3S+(1859/4104)*k4S-(11/40)*k5S));
    k6E = h*(b*(1-ts)*(S(i)-(8/27)*k1S+2*k2S-(3544/2565)*k3S+(1859/4104)*k4S-(11/40)*k5S)*((I(i)-(8/27)*k1I+2*k2I-(3544/2565)*k3I+(1859/4104)*k4I-(11/40)*k5I)+td*(E(i)-(8/27)*k1E+2*k2E-(3544/2565)*k3E+(1859/4104)*k4E-(11/40)*k5E))-(d+mu)*(E(i)-(8/27)*k1E+2*k2E-(3544/2565)*k3E+(1859/4104)*k4E-(11/40)*k5E));
    k6I = h*(d*(E(i)-(8/27)*k1E+2*k2E-(3544/2565)*k3E+(1859/4104)*k4E-(11/40)*k5E)-(g+s+mu)*(I(i)-(8/27)*k1I+2*k2I-(3544/2565)*k3I+(1859/4104)*k4I-(11/40)*k5I));
    k6R = h*((g+tt*s)*(I(i)-(8/27)*k1I+2*k2I-(3544/2565)*k3I+(1859/4104)*k4I-(11/40)*k5I)-mu*(R(i)-(8/27)*k1R+2*k2R-(3544/2565)*k3R+(1859/4104)*k4R-(11/40)*k5R));

S(i+1) = S(i)+(16/135)*k1S+(6656/12825)*k3S+(28561/56430)*k4S-(9/50)*k5S+(2/55)*k6S; 
E(i+1) = E(i)+(16/135)*k1E+(6656/12825)*k3E+(28561/56430)*k4E-(9/50)*k5E+(2/55)*k6E;
I(i+1) = I(i)+(16/135)*k1I+(6656/12825)*k3I+(28561/56430)*k4I-(9/50)*k5I+(2/55)*k6I;
R(i+1) = R(i)+(16/135)*k1R+(6656/12825)*k3R+(28561/56430)*k4R-(9/50)*k5R+(2/55)*k6R;
end

%Metode ODE45
options=odeset('RelTol',1e-13,'AbsTol',1e-13);
f=@(t,x)[L-b*(1-ts)*x(1)*(x(3)+td*x(2))-mu*x(1);b*(1-ts)*x(1)*(x(3)+td*x(2))-(d+mu)*x(2);d*x(2)-(g+s+mu)*x(3);(g+tt*s)*x(3)-mu*x(4)];
[t,xsol]=ode45(f,(0:h:tend),[S(1),E(1),I(1),R(1)],options);
 
%Gambar grafik
figure(1)
plot(t,S,'r',t,E,'g',t,I,'b',t,R,'c','linewidth',1.5);
xlabel('Waktu(t)'), ylabel('Nilai S(t),E(t),I(t),R(t)');
grid on
title('Penerapan Metode RKF45 Orde 5 Untuk theta_1 = 0.01')
legend('Susceptible (S(t))','Exposed (E(t))','Infected (I(t))','Recovered (R(t))');
 
figure(2)
plot(t,S,'b',t,xsol(:,1),'r--','linewidth',2);
xlabel('Waktu(t)'), ylabel('Nilai S(t)');
grid on
title('Perbandingan Metode RKF45 Orde 5 dengan ODE45 pada kelompok Susceptible')
legend('RKF45 orde 5','ODE45');
 
figure(3)
plot(t,E,'b',t,xsol(:,2),'r--','linewidth',2);
xlabel('Waktu(t)'), ylabel('Nilai E(t)');
grid on
title('Perbandingan Metode RKF45 Orde 5 dengan ODE45 pada kelompok Exposed')
legend('RKF45 orde 5','ODE45');
 
figure(4)
plot(t,I,'b',t,xsol(:,3),'r--','linewidth',2);
xlabel('Waktu(t)'), ylabel('Nilai I(t)');
grid on
title('Perbandingan Metode RKF45 Orde 5 dengan ODE45 pada kelompok Infected')
legend('RKF45 orde 5','ODE45');

figure(5)
plot(t,R,'b',t,xsol(:,4),'r--','linewidth',2);
xlabel('Waktu(t)'), ylabel('Nilai R(t)');
grid on
title('Perbandingan Metode RKF45 Orde 5 dengan ODE45 pada kelompok Recovered')
legend('RKF45 orde 5','ODE45');

%Error Mutlak
for i=1:N
    EMS(i)=abs(xsol((i),1)-S(i));
    EME(i)=abs(xsol((i),2)-E(i));
    EMI(i)=abs(xsol((i),3)-I(i));
    EMR(i)=abs(xsol((i),4)-R(i));
end

%Error Relatif
for i=1:N
    ES(i)=EMS(i)/abs(xsol((i),1));
    EE(i)=EME(i)/abs(xsol((i),2));
    EI(i)=EMI(i)/abs(xsol((i),3));
    ER(i)=EMR(i)/abs(xsol((i),4));
end

%Rerata Error Mutlak
rerataEMS=sum(EMS)/N
rerataEME=sum(EME)/N
rerataEMI=sum(EMI)/N
rerataEMR=sum(EMR)/N

%Rerata Error Relatif
rerataES=sum(ES)/N
rerataEE=sum(EE)/N
rerataEI=sum(EI)/N
rerataER=sum(ER)/N

bilanganreproduksi = (L*b*(1-ts)*(d+td*(g+s+mu)))/(mu*(g+s+mu)*(d+mu))