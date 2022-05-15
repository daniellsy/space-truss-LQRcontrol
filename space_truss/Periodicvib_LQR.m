clear;clc;
load para.mat;
mode=6;%截断模态数
num=size(gVout,2);%压电数目
v=v(:,1:mode);
M=v'*gM*v;
len=length(gM);
for i=1:mode
    v(:,i)=v(:,i)/sqrt(M(i,i));
end
M0=v'*gM*v;
K0=v'*gK*v;
w=zeros(mode,1);
kesi=zeros(mode,1);
for i=1:1:mode%set each mode damping
    kesi(i)=0.005+0.005*(i-1);
end
% kesi(2)=0.005;
for i=1:1:mode
    w(i)=sqrt(K0(i,i))/(2*pi);
end
C0=eye(mode);
for i=1:1:mode
    C0(i,i)=2*kesi(i)*w(i);
end
%增加阻尼项
Mn=eye(mode);Kn=zeros(mode);
for i=1:mode
    Kn(i,i)=K0(i,i);
end
Vout0=v'*gVout;
A=[zeros(mode) eye(mode);-Kn -C0];
B=[zeros(mode,num);v'*gVout];
% C=[gSout*v zeros(1,mode)];
Q=5000*eye(2*mode);
R=0.01*eye(num);
G1=lqr(A,B,Q,R);
A1=A-B*G1;
vp=84;%给定一点
vt=83;%观测点
% v=1e-2;
% x(size+(vp-1)*3+1)=v;
f=zeros(len,1);x=zeros(mode*2,1);x2=zeros(mode*2,1);
i=1;h=0.001;
for t=0:h:10
f1=50*sin(2*pi*11*t);
f((vp-1)*3+3)=f1;
F=[zeros(mode,1);v'*f];
K1=A*x+F;
K2=A*(x+h*K1/2)+F;
K3=A*(x+h*K2/2)+F;
K4=A*(x+h*K3)+F;
K12=A1*x2+F; 
K22=A1*(x2+h*K12/2)+F;
K32=A1*(x2+h*K22/2)+F;
K42=A1*(x2+h*K32)+F;
x=x+h/6*(K1+2*K2+2*K3+K4);
x2=x2+h/6*(K12+2*K22+2*K32+K42);
wy(i)=v((vt-1)*3+2,:)*x(1:mode,:);
wy2(i)=v((vt-1)*3+2,:)*x2(1:mode,:);
wz(i)=v((vt-1)*3+3,:)*x(1:mode,:);
wz2(i)=v((vt-1)*3+3,:)*x2(1:mode,:);
% Vs(i)=C*x;
u(:,i)=-G1*x2;
tc(i)=t;
i=i+1;
end
y=fft(wy);
ft=(0:length(y)-1)*1000/length(y);
plot(tc,wy,'--','LineWidth',0.6,'color','k');
hold on;
plot(tc,wy2,'LineWidth',1.4,'color','r');
xlabel('t/s');ylabel('w/m');
title('VibyMag');
legend('uncontrolled state','controlled state');
% figure(2)
% plot(tc,Vs,'LineWidth',1);
figure(2)
plot(tc,wz,'--','LineWidth',0.6,'color','k');
hold on;
plot(tc,wz2,'LineWidth',1.4,'color','r');
xlabel('t/s');ylabel('w/m');
title('VibzMag');
legend('uncontrolled state','controlled state');
figure(3)
plot(tc,u,'LineWidth',1);
xlabel('t/s');ylabel('u/V');
title('pzt-V');
legend('pzt-1','pzt-2','pzt-3','pzt-4');
figure(4)
plot(ft,abs(y));