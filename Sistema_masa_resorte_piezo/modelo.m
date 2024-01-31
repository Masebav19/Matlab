%% Respuesta en frecuencia
clear
clc
m1=157.6;
m2=54;
c1=42.98e3;
c2=40.11e3;
c0=1.89e-8;
k1=100.3e6;
k2=322.3e6;
alpha=1.52e-3;
R=30455.3;

Am=140e3;%amplitud de la entrada armónica

s=tf('s');

GV=(-alpha*R*s*(m1*s^2+k1+c1*s))/((m1*s^2+k1+c1*s)*((k2+c2*s+m2*s^2 )*(1+R*c0*s)+alpha^2*R*s)+...
    ((k2+c2*s)*(1+R*c0*s)+alpha^2*R*s)*m2*s^2);

GP=(-alpha^2*R*s^2*(m1*s^2+k1+c1*s)^2)/(2*((m1*s^2+k1+c1*s)*((k2+c2*s+m2*s^2 )*(1+R*c0*s)+alpha^2*R*s)+...
    ((k2+c2*s)*(1+R*c0*s)+alpha^2*R*s)*m2*s^2)^2);


A=[0 1 0 0 0;-(k1+k2)/m1 -(c1+c2)/m1 k2/m1 c2/m1 -alpha/m1;...
    0 0 0 1 0; k2/m2 c2/m2 -k2/m2 -c2/m2 alpha/m2;0 -alpha/c0 0 alpha/c0 -1/(R*c0)];
B=[0;0;0;1/m2;0];
C=[0 0 0 0 1];
D=0;

sys_V=ss(A,B,C,D);
Fs=linspace(0,800,100);
Vrms = zeros(size(Fs));
for i=1:length(Fs)
t=linspace(0,20,200000);
u=Am*sin(2*pi*Fs(i)*t);
y=lsim(sys_V,u,t);
Vrms(i)=sqrt(trapz(t,y.^2)/t(end));
end
figure
plot(Fs,Vrms)
hold on
Fs=linspace(0,800,100);
Vrms = zeros(size(Fs));
for i=1:length(Fs)
t=linspace(0,20,200000);
u=Am*sin(2*pi*Fs(i)*t);
y=lsim(GV,u,t);
Vrms(i)=sqrt(trapz(t,y.^2)/t(end));
end
plot(Fs,Vrms,'or')
title('Gráfica de respuesta del voltaje en función de la frecuencia');
xlabel('Frecuencia [Hz]');
ylabel('Voltaje [V]');
legend('Respuesta en el dominio del tiempo','Respuesta en el dominio de la frecuencia');
hold off

Fs=linspace(0,800,100);
P = zeros(size(Fs));
Vrms = zeros(size(Fs));
for i=1:length(Fs)
t=linspace(0,20,200000);
u=Am*sin(2*pi*Fs(i)*t);
y=lsim(sys_V,u,t);
Vrms(i)=sqrt(trapz(t,y.^2)/t(end));
P(i)=0.5*Vrms(i)^2/(Am*R/sqrt(2));
end
figure
plot(Fs,1000.*P)
hold on
Fs=linspace(0,800,100);
P = zeros(size(Fs));
for i=1:length(Fs)
t=linspace(0,20,200000);
u=Am*sin(2*pi*Fs(i)*t);
y=lsim(GP,u,t);
P(i)=sqrt(trapz(t,y.^2)/t(end));
end
plot(Fs,1000.*P,'or')
title('Gráfica de respuesta de la potencia en función de la frecuencia');
xlabel('Frecuencia [Hz]');
ylabel('Potencia [mW]');
legend('Respuesta en el dominio del tiempo','Respuesta en el dominio de la frecuencia');
hold off

%% variacion de parámetros del piezoeléctrico
% Capacitancia del piezo

m1=157.6;
m2=54;
c1=42.98e3;
c2=40.11e3;
c0=linspace(1.89e-8,1.89e-5,200);
k1=100.3e6;
k2=322.3e6;
alpha=1.52e-3;
R=30455.3;

Am=140e3;

s=tf('s');
Vrms = zeros(size(c0));
for i=1:length(c0)
    GV=(-alpha*R*s*(m1*s^2+k1+c1*s))/((m1*s^2+k1+c1*s)*((k2+c2*s+m2*s^2 )*(1+R*c0(i)*s)+alpha^2*R*s)+...
    ((k2+c2*s)*(1+R*c0(i)*s)+alpha^2*R*s)*m2*s^2);
t=linspace(0,20,200000);
u=Am*sin(2*pi*1.45*t);
y=lsim(GV,u,t);
Vrms(i)=sqrt(trapz(t,y.^2)/t(end));
end
figure
plot(c0,Vrms)
title('Respuesta del voltaje ante variación de parámetros');
xlabel('Capacitancia del piezoeléctrico [F]');
ylabel('Voltaje [V]');

% Factor de fuerza alpha
c0=1.89e-8;
alpha=linspace(1.52e-3,1.52e-1,200);

Am=10;

s=tf('s');
Vrms = zeros(size(alpha));
for i=1:length(alpha)
    GV=(-alpha(i)*R*s*(m1*s^2+k1+c1*s))/((m1*s^2+k1+c1*s)*((k2+c2*s+m2*s^2 )*(1+R*c0*s)+alpha(i)^2*R*s)+...
    ((k2+c2*s)*(1+R*c0*s)+alpha(i)^2*R*s)*m2*s^2);
t=linspace(0,20,200000);
u=Am*sin(2*pi*1.45*t);
y=lsim(GV,u,t);
Vrms(i)=sqrt(trapz(t,y.^2)/t(end));
end
figure
plot(alpha,Vrms)
title('Respuesta del voltaje ante variación de parámetros');
xlabel('Factor de fuerza [N/Volt]');
ylabel('Voltaje [V]');

% Resistencia

alpha=1.52e-3;
R=linspace(10e3,100e3,400);

Am=10;

s=tf('s');
Vrms = zeros(size(R));
for i=1:length(R)
    GV=(-alpha*R(i)*s*(m1*s^2+k1+c1*s))/((m1*s^2+k1+c1*s)*((k2+c2*s+m2*s^2 )*(1+R(i)*c0*s)+alpha^2*R(i)*s)+...
    ((k2+c2*s)*(1+R(i)*c0*s)+alpha^2*R(i)*s)*m2*s^2);
t=linspace(0,20,200000);
u=Am*sin(2*pi*1.45*t);
y=lsim(GV,u,t);
Vrms(i)=sqrt(trapz(t,y.^2)/t(end));
end
figure
plot(R,Vrms)
title('Respuesta del voltaje ante variación de parámetros');
xlabel('Resistencia [ohms]');
ylabel('Voltaje [V]');

%% tasa de masas fig 6.2

m1=linspace(0.1*157.6,500*157.6,100);
m2=54;
c1=42.98e3;
c2=40.11e3;
c0=1.89e-8;
k1=100.3e6;
k2=322.3e6;
alpha=1.52e-3;
R=30455.3;

Am=140e3;%amplitud de la entrada armónica

s=tf('s');

Fs=80;
P = zeros(size(m1));
Pin = zeros(size(m1));
for i=1:length(m1)
    
    A=[0 1 0 0 0;-(k1+k2)/m1(i) -(c1+c2)/m1(i) k2/m1(i) c2/m1(i) -alpha/m1(i);...
    0 0 0 1 0; k2/m2 c2/m2 -k2/m2 -c2/m2 alpha/m2;0 -alpha/c0 0 alpha/c0 -1/(R*c0)];
    B=[0;0;0;1/m2;0];
    C=[0 1 0 0 0;0 0 0 1 0];
    D=0;

    sys_V=ss(A,B,C,D);

    t=linspace(0,20,200000);
    u=Am*sin(2*pi*Fs*t);
    up=2*pi*Fs*Am*cos(2*pi*Fs*t);
    Y=lsim(sys_V,u,t);
    y1 = Y(:,1);
    dy1 = y1(2:end)-y1(1:end-1);
    y2 = Y(:,2);
    dy2 = y2(2:end)-y2(1:end-1);
    dt =  t(2:end)-t(1:end-1);
    Pi=0.5*(-m1(i).*(dy1'/dt).*up(1:end-1)')+0.5*(m2*(dy2'/dt).*up(1,end-1)');
    Pin(i)=sqrt(trapz(t(1,end-1),Pi.^2)/t(end-1));

    GP=(-alpha^2*R*s^2*(m1(i)*s^2+k1+c1*s)^2)/(2*((m1(i)*s^2+k1+c1*s)*((k2+c2*s+m2*s^2 )*(1+R*c0*s)+alpha^2*R*s)+...
    ((k2+c2*s)*(1+R*c0*s)+alpha^2*R*s)*m2*s^2)^2);
    t=linspace(0,20,200000);
    u=Am*sin(2*pi*Fs*t);
    y=lsim(GP,u,t);
    P(i)=sqrt(trapz(t,y.^2)/t(end));
end
efic = P./Pin;
figure
subplot(2,2,1)
plot(m1./m2,1000.*P,'b')
title('Respuesta de la potencia')
xlabel('Tasa de cambio de masas m1/m2')
ylabel('Potencia [mW]')
legend('Variación de m1')

subplot(2,2,2)
plot(m1./m2,efic,'b')
title('Respuesta de la eficiencia')
xlabel('Tasa de cambio de masas m1/m2')
ylabel('Eficiencia')
legend('Variación de m1')

m2= linspace(0.1*54,5*54,100);
m1=m1(end);

Fs=80;
P = zeros(size(m2));
Pin = zeros(size(m2));
for i=1:length(m2)
    
    A=[0 1 0 0 0;-(k1+k2)/m1 -(c1+c2)/m1 k2/m1 c2/m1 -alpha/m1;...
    0 0 0 1 0; k2/m2(i) c2/m2(i) -k2/m2(i) -c2/m2(i) alpha/m2(i);0 -alpha/c0 0 alpha/c0 -1/(R*c0)];
    B=[0;0;0;1/m2(i);0];
    C=[0 1 0 0 0;0 0 0 1 0];
    D=0;

    sys_V=ss(A,B,C,D);

    t=linspace(0,20,200000);
    u=Am*sin(2*pi*Fs*t);
    up=2*pi*Fs*Am*cos(2*pi*Fs*t);
    Y=lsim(sys_V,u,t);
    y1 = Y(:,1);
    dy1 = y1(2:end)-y1(1:end-1);
    y2 = Y(:,2);
    dy2 = y2(2:end)-y2(1:end-1);
    dt =  t(2:end)-t(1:end-1);
    Pi=0.5*(-m1.*(dy1'/dt).*up(1:end-1)')+0.5*(m2(i)*(dy2'/dt).*up(1,end-1)');
    Pin(i)=sqrt(trapz(t(1,end-1),Pi.^2)/t(end-1));
    
    GP=(-alpha^2*R*s^2*(m1*s^2+k1+c1*s)^2)/(2*((m1*s^2+k1+c1*s)*((k2+c2*s+m2(i)*s^2 )*(1+R*c0*s)+alpha^2*R*s)+...
    ((k2+c2*s)*(1+R*c0*s)+alpha^2*R*s)*m2(i)*s^2)^2);
t=linspace(0,20,200000);
u=Am*sin(2*pi*Fs*t);
y=lsim(GP,u,t);
P(i)=sqrt(trapz(t,y.^2)/t(end));
end
efic = P./Pin;
subplot(2,2,3)
plot(m1./m2,1000.*P,'r')
title('Respuesta de la potencia')
xlabel('Tasa de cambio de masas m1/m2')
ylabel('Potencia [mW]')
legend('Variación de m2')

subplot(2,2,4)
plot(m1./m2,efic,'r')
title('Respuesta de la eficiencia')
xlabel('Tasa de cambio de masas m1/m2')
ylabel('Eficiencia')
legend('Variación de m2')

%% tasa de rigidez fig 6.3

m1=157.6;
m2=54;
c1=42.98e3;
c2=40.11e3;
c0=1.89e-8;
k1=linspace(0.1*100.3e6,50*100.3e6,100);
k2=322.3e6;
alpha=1.52e-3;
R=30455.3;

Am=140e3;%amplitud de la entrada armónica

s=tf('s');

Fs=80;
P = zeros(size(k1));
Pin = zeros(size(k1));
for i=1:length(k1)
    
    A=[0 1 0 0 0;-(k1(i)+k2)/m1 -(c1+c2)/m1 k2/m1 c2/m1 -alpha/m1;...
    0 0 0 1 0; k2/m2 c2/m2 -k2/m2 -c2/m2 alpha/m2;0 -alpha/c0 0 alpha/c0 -1/(R*c0)];
    B=[0;0;0;1/m2;0];
    C=[0 1 0 0 0;0 0 0 1 0];
    D=0;

    sys_V=ss(A,B,C,D);

    t=linspace(0,10,200000);
    u=Am*sin(2*pi*Fs*t);
    up=2*pi*Fs*Am*cos(2*pi*Fs*t);
    Y=lsim(sys_V,u,t);
    y1 = Y(:,1);
    dy1 = y1(2:end)-y1(1:end-1);
    y2 = Y(:,2);
    dy2 = y2(2:end)-y2(1:end-1);
    dt =  t(2:end)-t(1:end-1);
    Pi=0.5*(-m1.*(dy1'/dt).*up(1:end-1)')+0.5*(m2*(dy2'/dt).*up(1,end-1)');
    Pin(i)=sqrt(trapz(t(1,end-1),Pi.^2)/t(end-1));

    GP=(-alpha^2*R*s^2*(m1*s^2+k1(i)+c1*s)^2)/(2*((m1*s^2+k1(i)+c1*s)*((k2+c2*s+m2*s^2 )*(1+R*c0*s)+alpha^2*R*s)+...
    ((k2+c2*s)*(1+R*c0*s)+alpha^2*R*s)*m2*s^2)^2);
    t=linspace(0,10,200000);
    u=Am*sin(2*pi*Fs*t);
    y=lsim(GP,u,t);
    P(i)=sqrt(trapz(t,y.^2)/t(end));
end
efic = P./Pin;
figure
subplot(2,2,1)
plot(k1./k2,1000.*P,'b')
title('Respuesta de la potencia')
xlabel('Tasa de cambio de rigidez k1/k2')
ylabel('Potencia [mW]')
legend('Variación de k1')

subplot(2,2,2)
plot(k1./k2,efic,'b')
title('Respuesta de la eficiencia')
xlabel('Tasa de cambio de rigidez k1/k2')
ylabel('Eficiencia')
legend('Variación de k1')


k1=k1(end);
k2=linspace(0.1*322.3e6,50*322.3e6,100);

Fs=80;
P = zeros(size(k2));
Pin = zeros(size(k2));
for i=1:length(k2)
    
    A=[0 1 0 0 0;-(k1+k2(i))/m1 -(c1+c2)/m1 k2(i)/m1 c2/m1 -alpha/m1;...
    0 0 0 1 0; k2(i)/m2 c2/m2 -k2(i)/m2 -c2/m2 alpha/m2;0 -alpha/c0 0 alpha/c0 -1/(R*c0)];
    B=[0;0;0;1/m2;0];
    C=[0 1 0 0 0;0 0 0 1 0];
    D=0;

    sys_V=ss(A,B,C,D);

    t=linspace(0,10,200000);
    u=Am*sin(2*pi*Fs*t);
    up=2*pi*Fs*Am*cos(2*pi*Fs*t);
    Y=lsim(sys_V,u,t);
    y1 = Y(:,1);
    dy1 = y1(2:end)-y1(1:end-1);
    y2 = Y(:,2);
    dy2 = y2(2:end)-y2(1:end-1);
    dt =  t(2:end)-t(1:end-1);
    Pi=0.5*(-m1.*(dy1'/dt).*up(1:end-1)')+0.5*(m2*(dy2'/dt).*up(1,end-1)');
    Pin(i)=sqrt(trapz(t(1,end-1),Pi.^2)/t(end-1));
    
    GP=(-alpha^2*R*s^2*(m1*s^2+k1+c1*s)^2)/(2*((m1*s^2+k1+c1*s)*((k2(i)+c2*s+m2*s^2 )*(1+R*c0*s)+alpha^2*R*s)+...
    ((k2(i)+c2*s)*(1+R*c0*s)+alpha^2*R*s)*m2*s^2)^2);
t=linspace(0,10,200000);
u=Am*sin(2*pi*Fs*t);
y=lsim(GP,u,t);
P(i)=sqrt(trapz(t,y.^2)/t(end));
end
efic = P./Pin;
subplot(2,2,3)
plot(k1./k2,1000.*P,'r')
title('Respuesta de la potencia')
xlabel('Tasa de cambio de regidez k1/k2')
ylabel('Potencia [mW]')
legend('Variación de k2')

subplot(2,2,4)
plot(k1./k2,efic,'r')
title('Respuesta de la eficiencia')
xlabel('Tasa de cambio de rigidez k1/k2')
ylabel('Eficiencia')
legend('Variación de k2')

%% tasa de desplazamiento fig 6.8

m1=157.6;
m2=54;
c1=42.98e3;
c2=40.11e3;
c0=1.89e-8;
k1=100.3e6;
k2=322.3e6;
alpha=1.52e-3;
R=30455.3;

Am=140e3;%amplitud de la entrada armónica

Fs=linspace(0.1,800,100);
x1=zeros(size(Fs));
x2=zeros(size(Fs));
Y1=zeros(size(Fs));
Y2=zeros(size(Fs));
for i=1:length(Fs)
    
    A=[0 1 0 0 0;-(k1+k2)/m1 -(c1+c2)/m1 k2/m1 c2/m1 -alpha/m1;...
    0 0 0 1 0; k2/m2 c2/m2 -k2/m2 -c2/m2 alpha/m2;0 -alpha/c0 0 alpha/c0 -1/(R*c0)];
    B=[0;0;0;1/m2;0];
    C=[1 0 0 0 0;0 0 1 0 0];
    D=0;

    sys_V=ss(A,B,C,D);

    t=linspace(0,10,200000);
    u=Am*sin(2*pi*Fs(i)*t);
    Y=lsim(sys_V,u,t);
    Y1(i)=sqrt(trapz(t,Y(:,1).^2)/t(end));
    Y2(i)=sqrt(trapz(t,Y(:,2).^2)/t(end));
    x1(i) = sqrt(2)*Y1(i)./Am;
    x2(i) = sqrt(2)*Y2(i)./Am;
end
figure
plot(Fs,x1,'b');
hold on
plot(Fs,x2,'--r');
title('Respuesta de tasa de desplazamiento');
xlabel('Frecuencia [Hz]')
ylabel('Tasa de desplzamiento');
legend('X1/Y','X2/Y');
hold off
figure
plot(Fs,Y1./Y2);
title('Respuesta de tasa de desplazamiento');
xlabel('Frecuencia [Hz]')
ylabel('Tasa de desplzamiento X1/X2');

%% Variación de resistencia fig 6.11

m1=157.6;
m2=54;
c1=42.98e3;
c2=40.11e3;
c0=1.89e-8;
k1=100.3e6;
k2=322.3e6;
alpha=1.52e-3;
R=linspace(0.1*30455.3,30*30455.3,100);

Am=140e3;%amplitud de la entrada armónica

s=tf('s');

Fs=444;
P = zeros(size(R));
for i=1:length(R)
    GP=(-alpha^2*R(i)*s^2*(m1*s^2+k1+c1*s)^2)/(2*((m1*s^2+k1+c1*s)*((k2+c2*s+m2*s^2 )*(1+R(i)*c0*s)+alpha^2*R(i)*s)+...
    ((k2+c2*s)*(1+R(i)*c0*s)+alpha^2*R(i)*s)*m2*s^2)^2);
    t=linspace(0,20,200000);
    u=Am*sin(2*pi*Fs*t);
    y=lsim(GP,u,t);
    P(i)=sqrt(trapz(t,y.^2)/t(end));
end
figure
plot(R,1000.*P,'b')
title('Respuesta de la potencia')
xlabel('Resistencia [ohmios]')
ylabel('Potencia [mW]')

%% variación de c1 fig 6.12
m1=157.6;
m2=54;
c1=linspace(0.1*42.98e3,1e2*42.98e3,100);
c2=40.11e3;
c0=1.89e-8;
k1=100.3e6;
k2=322.3e6;
alpha=1.52e-3;
R=30455.3;

Am=140e3;%amplitud de la entrada armónica

s=tf('s');

Fs=150;
P = zeros(size(c1));
for i=1:length(c1)
    GP=(-alpha^2*R*s^2*(m1*s^2+k1+c1(i)*s)^2)/(2*((m1*s^2+k1+c1(i)*s)*((k2+c2*s+m2*s^2 )*(1+R*c0*s)+alpha^2*R*s)+...
    ((k2+c2*s)*(1+R*c0*s)+alpha^2*R*s)*m2*s^2)^2);
    t=linspace(0,20,200000);
    u=Am*sin(2*pi*Fs*t);
    y=lsim(GP,u,t);
    P(i)=sqrt(trapz(t,y.^2)/t(end));
end
figure
plot(c1./42.98e3,1000.*P,'b')
title('Respuesta de la potencia')
xlabel('coeficiente de amortiguamiento [x42.98e3 N*s/m]')
ylabel('Potencia [mW]')

%% coeficiente de amortiguamiento de suspensión fig 6.13
m1=157.6;
m2=54;
c1=42.98e3;
c2=linspace(0.1*40.11e3,1e2*40.11e3,100);
c0=1.89e-8;
k1=100.3e6;
k2=322.3e6;
alpha=1.52e-3;
R=30455.3;

Am=140e3;%amplitud de la entrada armónica

s=tf('s');

Fs=150;
P = zeros(size(c2));
for i=1:length(c2)
    GP=(-alpha^2*R*s^2*(m1*s^2+k1+c1*s)^2)/(2*((m1*s^2+k1+c1*s)*((k2+c2(i)*s+m2*s^2 )*(1+R*c0*s)+alpha^2*R*s)+...
    ((k2+c2(i)*s)*(1+R*c0*s)+alpha^2*R*s)*m2*s^2)^2);
    t=linspace(0,20,200000);
    u=Am*sin(2*pi*Fs*t);
    y=lsim(GP,u,t);
    P(i)=sqrt(trapz(t,y.^2)/t(end));
end
figure
plot(c2./40.11e3,1000.*P,'b')
title('Respuesta de la potencia')
xlabel('coeficiente de amortiguamiento [x40.11e3 N*s/m]')
ylabel('Potencia [mW]')

%% factor de fuerza fig 6.14

m1=157.6;
m2=54;
c1=42.98e3;
c2=40.11e3;
c0=1.89e-8;
k1=100.3e6;
k2=322.3e6;
alpha=linspace(0.1*1.52e-3,2e4*1.52e-3,100);
R=30455.3;

Am=140e3;%amplitud de la entrada armónica

s=tf('s');

Fs=150;
P = zeros(size(alpha));
for i=1:length(alpha)
    GP=(-alpha(i)^2*R*s^2*(m1*s^2+k1+c1*s)^2)/(2*((m1*s^2+k1+c1*s)*((k2+c2*s+m2*s^2 )*(1+R*c0*s)+alpha(i)^2*R*s)+...
    ((k2+c2*s)*(1+R*c0*s)+alpha(i)^2*R*s)*m2*s^2)^2);
    t=linspace(0,10,2000000);
    u=Am*sin(2*pi*Fs*t);
    y=lsim(GP,u,t);
    P(i)=sqrt(trapz(t,y.^2)/t(end));
end
figure
plot(alpha./1.52e-3,1000.*P,'b')
title('Respuesta de la potencia')
xlabel('Factor de fuerza [x1.52e-3 N/volt]')
ylabel('Potencia [mW]')

%% variacion de masa 1 fig 6.15

m1=[120 157.6 200 270];
m2=54;
c1=42.98e3;
c2=40.11e3;
c0=1.89e-8;
k1=100.3e6;
k2=322.3e6;
alpha=1.52e-3;
R=30455.3;

Am=140e3;%amplitud de la entrada armónica

s=tf('s');

Fs=linspace(0,800,150);
P = zeros(size(Fs));
figure
hold on
for j=1:length(m1)
for i=1:length(Fs)
    GP=(-alpha^2*R*s^2*(m1(j)*s^2+k1+c1*s)^2)/(2*((m1(j)*s^2+k1+c1*s)*((k2+c2*s+m2*s^2 )*(1+R*c0*s)+alpha^2*R*s)+...
    ((k2+c2*s)*(1+R*c0*s)+alpha^2*R*s)*m2*s^2)^2);
    t=linspace(0,10,2000000);
    u=Am*sin(2*pi*Fs(i)*t);
    y=lsim(GP,u,t);
    P(i)=sqrt(trapz(t,y.^2)/t(end));
end
msg=sprintf("m1: %.2f",m1(j));
plot(Fs,1000.*P,'DisplayName',msg)
end
title('Respuesta de la potencia')
xlabel('Frecuencia [Hz]')
ylabel('Potencia [mW]')
legend
hold off

%% Variación masa 2 fig 6.16

m1=157.6;
m2=[54 60];
c1=42.98e3;
c2=40.11e3;
c0=1.89e-8;
k1=100.3e6;
k2=322.3e6;
alpha=1.52e-3;
R=30455.3;

Am=140e3;%amplitud de la entrada armónica

s=tf('s');

Fs=linspace(0,800,150);
P = zeros(size(Fs));
figure
hold on
for j=1:length(m2)
for i=1:length(Fs)
    GP=(-alpha^2*R*s^2*(m1*s^2+k1+c1*s)^2)/(2*((m1*s^2+k1+c1*s)*((k2+c2*s+m2(j)*s^2 )*(1+R*c0*s)+alpha^2*R*s)+...
    ((k2+c2*s)*(1+R*c0*s)+alpha^2*R*s)*m2(j)*s^2)^2);
    t=linspace(0,10,2000000);
    u=Am*sin(2*pi*Fs(i)*t);
    y=lsim(GP,u,t);
    P(i)=sqrt(trapz(t,y.^2)/t(end));
end
msg=sprintf("m2: %.2f",m2(j));
plot(Fs,1000.*P,'DisplayName',msg)
end
title('Respuesta de la potencia')
xlabel('Frecuencia [Hz]')
ylabel('Potencia [mW]')
legend
hold off
%% Variación de k1 fig 6.17

m1=157.6;
m2=54;
c1=42.98e3;
c2=40.11e3;
c0=1.89e-8;
k1=[22.89e6 83.3e6 100.3e6 257.9e6];
k2=322.3e6;
alpha=1.52e-3;
R=30455.3;

Am=140e3;%amplitud de la entrada armónica

s=tf('s');

Fs=linspace(0,800,150);
P = zeros(size(Fs));
figure
hold on
for j=1:length(k1)
for i=1:length(Fs)
    GP=(-alpha^2*R*s^2*(m1*s^2+k1(j)+c1*s)^2)/(2*((m1*s^2+k1(j)+c1*s)*((k2+c2*s+m2*s^2 )*(1+R*c0*s)+alpha^2*R*s)+...
    ((k2+c2*s)*(1+R*c0*s)+alpha^2*R*s)*m2*s^2)^2);
    t=linspace(0,10,2000000);
    u=Am*sin(2*pi*Fs(i)*t);
    y=lsim(GP,u,t);
    P(i)=sqrt(trapz(t,y.^2)/t(end));
end
msg=sprintf("k1: %.2fe6",k1(j)/1e6);
plot(Fs,1000.*P,'DisplayName',msg)
end
title('Respuesta de la potencia')
xlabel('Frecuencia [Hz]')
ylabel('Potencia [mW]')
legend
hold off

%% Variación de k2 fig 6.18

m1=157.6;
m2=54;
c1=42.98e3;
c2=40.11e3;
c0=1.89e-8;
k1=100.3e6;
k2=[192.6e6 322.3e6 401.1e6 583.3e6];
alpha=1.52e-3;
R=30455.3;

Am=140e3;%amplitud de la entrada armónica

s=tf('s');

Fs=linspace(0,800,150);
P = zeros(size(Fs));
figure
hold on
for j=1:length(k2)
for i=1:length(Fs)
    GP=(-alpha^2*R*s^2*(m1*s^2+k1+c1*s)^2)/(2*((m1*s^2+k1+c1*s)*((k2(j)+c2*s+m2*s^2 )*(1+R*c0*s)+alpha^2*R*s)+...
    ((k2(j)+c2*s)*(1+R*c0*s)+alpha^2*R*s)*m2*s^2)^2);
    t=linspace(0,10,2000000);
    u=Am*sin(2*pi*Fs(i)*t);
    y=lsim(GP,u,t);
    P(i)=sqrt(trapz(t,y.^2)/t(end));
end
msg=sprintf("k2: %.2fe6",k2(j)/1e6);
plot(Fs,1000.*P,'DisplayName',msg)
end
title('Respuesta de la potencia')
xlabel('Frecuencia [Hz]')
ylabel('Potencia [mW]')
legend
hold off

%% variación de c1 fig 6.20

m1=157.6;
m2=54;
c1=[42.98e3 80.96e3 155.32e3 206e3];
c2=40.11e3;
c0=1.89e-8;
k1=100.3e6;
k2=322.3e6;
alpha=1.52e-3;
R=30455.3;

Am=140e3;%amplitud de la entrada armónica

s=tf('s');

Fs=linspace(0,800,150);
P = zeros(size(Fs));
figure
hold on
for j=1:length(c1)
for i=1:length(Fs)
    GP=(-alpha^2*R*s^2*(m1*s^2+k1+c1(j)*s)^2)/(2*((m1*s^2+k1+c1(j)*s)*((k2+c2*s+m2*s^2 )*(1+R*c0*s)+alpha^2*R*s)+...
    ((k2+c2*s)*(1+R*c0*s)+alpha^2*R*s)*m2*s^2)^2);
    t=linspace(0,10,2000000);
    u=Am*sin(2*pi*Fs(i)*t);
    y=lsim(GP,u,t);
    P(i)=sqrt(trapz(t,y.^2)/t(end));
end
msg=sprintf("c1: %.2fe3",c1(j)/1e3);
plot(Fs,1000.*P,'DisplayName',msg)
end
title('Respuesta de la potencia')
xlabel('Frecuencia [Hz]')
ylabel('Potencia [mW]')
legend
hold off

%% variación de c1 fig 6.21

m1=157.6;
m2=54;
c1=42.98e3;
c2=[40.11e3 55.14e3 70.21e3 98.26e3];
c0=1.89e-8;
k1=100.3e6;
k2=322.3e6;
alpha=1.52e-3;
R=30455.3;

Am=140e3;%amplitud de la entrada armónica

s=tf('s');

Fs=linspace(0,800,150);
P = zeros(size(Fs));
figure
hold on
for j=1:length(c2)
for i=1:length(Fs)
    GP=(-alpha^2*R*s^2*(m1*s^2+k1+c1*s)^2)/(2*((m1*s^2+k1+c1*s)*((k2+c2(j)*s+m2*s^2 )*(1+R*c0*s)+alpha^2*R*s)+...
    ((k2+c2(j)*s)*(1+R*c0*s)+alpha^2*R*s)*m2*s^2)^2);
    t=linspace(0,10,2000000);
    u=Am*sin(2*pi*Fs(i)*t);
    y=lsim(GP,u,t);
    P(i)=sqrt(trapz(t,y.^2)/t(end));
end
msg=sprintf("c2: %.2fe3",c2(j)/1e3);
plot(Fs,1000.*P,'DisplayName',msg)
end
title('Respuesta de la potencia')
xlabel('Frecuencia [Hz]')
ylabel('Potencia [mW]')
legend
hold off

%% tasa de amortiguamiento fig 6.22

m1=157.6;
m2=54;
c1=linspace(0.1*42.98e3,50*42.98e3,100);
c2=40.11e3;
c0=1.89e-8;
k1=100.3e6;
k2=322.3e6;
alpha=1.52e-3;
R=30455.3;

Am=140e3;%amplitud de la entrada armónica

s=tf('s');

Fs=444;
P = zeros(size(c1));
Pin = zeros(size(c1));
for i=1:length(c1)
    
    A=[0 1 0 0 0;-(k1+k2)/m1 -(c1(i)+c2)/m1 k2/m1 c2/m1 -alpha/m1;...
    0 0 0 1 0; k2/m2 c2/m2 -k2/m2 -c2/m2 alpha/m2;0 -alpha/c0 0 alpha/c0 -1/(R*c0)];
    B=[0;0;0;1/m2;0];
    C=[0 1 0 0 0;0 0 0 1 0];
    D=0;

    sys_V=ss(A,B,C,D);

    t=linspace(0,10,200000);
    u=Am*sin(2*pi*Fs*t);
    up=2*pi*Fs*Am*cos(2*pi*Fs*t);
    Y=lsim(sys_V,u,t);
    y1 = Y(:,1);
    dy1 = y1(2:end)-y1(1:end-1);
    y2 = Y(:,2);
    dy2 = y2(2:end)-y2(1:end-1);
    dt =  t(2:end)-t(1:end-1);
    Pi=0.5*(-m1.*(dy1'/dt).*up(1:end-1)')+0.5*(m2*(dy2'/dt).*up(1,end-1)');
    Pin(i)=sqrt(trapz(t(1,end-1),Pi.^2)/t(end-1));

    GP=(-alpha^2*R*s^2*(m1*s^2+k1+c1(i)*s)^2)/(2*((m1*s^2+k1+c1(i)*s)*((k2+c2*s+m2*s^2 )*(1+R*c0*s)+alpha^2*R*s)+...
    ((k2+c2*s)*(1+R*c0*s)+alpha^2*R*s)*m2*s^2)^2);
    t=linspace(0,10,200000);
    u=Am*sin(2*pi*Fs*t);
    y=lsim(GP,u,t);
    P(i)=sqrt(trapz(t,y.^2)/t(end));
end
efic = P./Pin;
figure
subplot(2,2,1)
plot(c1./c2,1000.*P,'b')
title('Respuesta de la potencia')
xlabel('Tasa deamortiguamiento c1/c2')
ylabel('Potencia [mW]')
legend('Variación de c1')

subplot(2,2,2)
plot(c1./c2,efic,'b')
title('Respuesta de la eficiencia')
xlabel('Tasa deamortiguamiento c1/c2')
ylabel('Eficiencia')
legend('Variación de c1')


c1=c1(end);
c2=linspace(0.1*40.11e3,50*40.11e3,100);

Fs=444;
P = zeros(size(c2));
Pin = zeros(size(c2));
for i=1:length(c2)
    
    A=[0 1 0 0 0;-(k1+k2)/m1 -(c1+c2(i))/m1 k2/m1 c2(i)/m1 -alpha/m1;...
    0 0 0 1 0; k2/m2 c2(i)/m2 -k2/m2 -c2(i)/m2 alpha/m2;0 -alpha/c0 0 alpha/c0 -1/(R*c0)];
    B=[0;0;0;1/m2;0];
    C=[0 1 0 0 0;0 0 0 1 0];
    D=0;

    sys_V=ss(A,B,C,D);

    t=linspace(0,10,200000);
    u=Am*sin(2*pi*Fs*t);
    up=2*pi*Fs*Am*cos(2*pi*Fs*t);
    Y=lsim(sys_V,u,t);
    y1 = Y(:,1);
    dy1 = y1(2:end)-y1(1:end-1);
    y2 = Y(:,2);
    dy2 = y2(2:end)-y2(1:end-1);
    dt =  t(2:end)-t(1:end-1);
    Pi=0.5*(-m1.*(dy1'/dt).*up(1:end-1)')+0.5*(m2*(dy2'/dt).*up(1,end-1)');
    Pin(i)=sqrt(trapz(t(1,end-1),Pi.^2)/t(end-1));
    
    GP=(-alpha^2*R*s^2*(m1*s^2+k1+c1*s)^2)/(2*((m1*s^2+k1+c1*s)*((k2+c2(i)*s+m2*s^2 )*(1+R*c0*s)+alpha^2*R*s)+...
    ((k2+c2(i)*s)*(1+R*c0*s)+alpha^2*R*s)*m2*s^2)^2);
t=linspace(0,10,200000);
u=Am*sin(2*pi*Fs*t);
y=lsim(GP,u,t);
P(i)=sqrt(trapz(t,y.^2)/t(end));
end
efic = P./Pin;
subplot(2,2,3)
plot(c1./c2,1000.*P,'r')
title('Respuesta de la potencia')
xlabel('Tasa deamortiguamiento c1/c2')
ylabel('Potencia [mW]')
legend('Variación de c2')

subplot(2,2,4)
plot(c1./c2,efic,'r')
title('Respuesta de la eficiencia')
xlabel('Tasa deamortiguamiento c1/c2')
ylabel('Eficiencia')
legend('Variación de c2')