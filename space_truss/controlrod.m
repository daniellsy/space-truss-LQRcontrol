function [M,Kn,Ku] = controlrod()
global gpMaterial
roup=gpMaterial(1);Ap=gpMaterial(2);Ep=gpMaterial(3);lp=gpMaterial(4);
np=gpMaterial(5);d33=gpMaterial(6);miup=gpMaterial(7);rou=gpMaterial(8);
A=gpMaterial(9);E=gpMaterial(10);l=gpMaterial(11);
m1=rou*A*l/6;
m3=rou*A*l/6;
m2=roup*Ap*lp/6;
k1=E*A/l;
k3=E*A/l;
k2=Ep*Ap/lp;
kuf=np*d33*Ap*Ep/lp;
Muu=[2*m1 m1 0 0;m1 2*m1+2*m2 m2 0;0 m2 2*m2+2*m3 m3;0 0 m3 2*m3];
Kuu=[k1 -k1 0 0;-k1 k1+k2 -k2 0;0 -k2 k2+k3 -k3;0 0 -k3 k3];
Kuf=kuf*[0;-1;1;0];
K=k1*k2+k2*k3+k1*k3;
Q=1/K*[K 0;k1*k2+k1*k3 k2*k3;k1*k2 k2*k3+k1*k3;0 K];
R=1/K*[0;k3*kuf;-k1*kuf;0];
Ms=Q'*Muu*Q;
Ks=Q'*Kuu*Q;
Kn=Ks(1,1);
Kus=Q'*(Kuu*R+Kuf);
L=l+lp+l;
Qv=1/L*[L 0;lp+l l;l l+lp;0 L];
M2s=Qv'*Muu*Qv;
M3s=Qv'*Muu*Qv; 
M=zeros(6);Ku=zeros(6,1);
M(1,1)=Ms(1,1);M(1,4)=Ms(1,2);M(4,1)=Ms(2,1);M(4,4)=Ms(2,2);
M(2,2)=M2s(1,1);M(2,5)=M2s(1,2);M(5,2)=M2s(2,1);M(5,5)=M2s(2,2);
M(3,3)=M3s(1,1);M(3,6)=M3s(1,2);M(6,3)=M3s(2,1);M(6,6)=M3s(2,2);
Ku(1)=Kus(1);Ku(4)=Kus(2);
end

