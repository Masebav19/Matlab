
clear
clc
m1=40;
m2=260;
c1=264.73;
c2=520;
c0=1.89e-8;
k1=130e3;
k2=26e3;
alpha=1.52e-3;
fn=1.45;
R=30455.3;

s=tf('s');

GV=(-m2*alpha*R*s*(k1+c1*s))/((m1*s^2+k1+c1*s)*((k2+c2*s+m2*s^2 )*(1+R*c0*s)+alpha^2*R*s)+...
    ((k2+c2*s)*(1+R*c0*s)+alpha^2*R*s)*m2*s^2);

Fs=linspace(0,14,200);
Vrms = zeros(size(Fs));
for i=1:length(Fs)
t=linspace(0,20,2000);
u=10*sin(2*pi*Fs(i)*t);
y=lsim(GV,u,t);
Vrms(i)=sqrt(trapz(t,y.^2)/t(end));
end
figure
plot(Fs,Vrms)