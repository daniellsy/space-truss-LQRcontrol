%    ����ƽ���ϵ������Ԫģ�����ݣ�
%        gNode ------- ��㶨��
%        gElement ---- ��Ԫ����
%        gMaterial --- ���϶��壬��������ģ�������Ľ���� 
%        gBC --------- Լ������ 
%        gK  --------- ����նȾ���
%        
clear;clc;tic;
global gNode gElement gMaterial gpMaterial gBC gK gM gVout cont
cont=1;pztn=4;%ѹ�������˵�Ԫ��
% ��������
%             ����ģ��   �����   �ܶ�  ���ɱ�
% gMaterial = [7.0e10   2.0e-4   2700   0.3       % ���� 1 (����
%              7.0e10   2.0e-4   2700   0.3       % ���� 2 (����
%              2.1e11   3.0e-4   7800   0.3];     % ���� 3���֣�
gMaterial = [7.27e10   6.283e-5   3100   0.3       % ��˲��� 1 (����
             7.27e10   6.283e-5   3100   0.3       % ���˲��� 2 (����
             7.27e10   6.283e-5   3100   0.3];     % б�˲��� 3���֣�
gpMaterial=[7600 7e-4 6.3e10 0.16 150 6.45e-10 1.5e-8 8000 6.283e-5 21e10 0.1];
%           roup    Ap     Ep    lp   np   d33     miup   rou    A       E    l
dx=0.3;%�쳤�����ϵ�һ�ڵĳ���
dy=0.4;%��ȷ�����һ�ڵĳ���
m=21;%����
mn=0;%��ͷ����
mtop=0;
gNode=zeros((m+1)*4,3);    
for i=1:m+1
      gNode((i-1)*4+1,1:3)=[dx*(i-1),0,0];  %�ڵ����꣨x,y,z��
      gNode((i-1)*4+2,1:3)=[dx*(i-1),dy,0];  %�ڵ����꣨x,y,z��
      gNode((i-1)*4+3,1:3)=[dx*(i-1),dy,dy];  %�ڵ����꣨x,y,z��
      gNode((i-1)*4+4,1:3)=[dx*(i-1),0,dy];  %�ڵ����꣨x,y,z��
end
 % ���嵥Ԫ
gnum=13;%����ÿ����Ԫ�м�����
gElement=zeros(m*gnum,3);
for i=1:m
        gElement((i-1)*gnum+1,1:3)=[(i-1)*4+1,(i-1)*4+5,1];%���Ϊ1 ����Ϊ2 б��Ϊ3
        gElement((i-1)*gnum+2,1:3)=[(i-1)*4+1,(i-1)*4+8,3]; 
        gElement((i-1)*gnum+3,1:3)=[(i-1)*4+2,(i-1)*4+6,1]; 
        gElement((i-1)*gnum+4,1:3)=[(i-1)*4+2,(i-1)*4+5,3]; 
        gElement((i-1)*gnum+5,1:3)=[(i-1)*4+3,(i-1)*4+7,1]; 
        gElement((i-1)*gnum+6,1:3)=[(i-1)*4+3,(i-1)*4+6,3]; 
        gElement((i-1)*gnum+7,1:3)=[(i-1)*4+4,(i-1)*4+8,1]; 
        gElement((i-1)*gnum+8,1:3)=[(i-1)*4+4,(i-1)*4+7,3];  
        gElement((i-1)*gnum+9,1:3)=[(i-1)*4+5,(i-1)*4+6,2]; 
        gElement((i-1)*gnum+10,1:3)=[(i-1)*4+5,(i-1)*4+8,2]; 
        gElement((i-1)*gnum+11,1:3)=[(i-1)*4+6,(i-1)*4+7,2]; 
        gElement((i-1)*gnum+12,1:3)=[(i-1)*4+7,(i-1)*4+8,2];  
        gElement((i-1)*gnum+13,1:3)=[(i-1)*4+6,(i-1)*4+8,3]; 
end
     % Լ������
     %     ����   ���ɶȺ�    Լ��ֵ 
     gBC = [ 1         1          0 
             1         2          0
             1         3          0
             2         1          0 
             2         2          0
             2         3          0
             3         1          0
             3         2          0
             3         3          0
             4         1          0
             4         2          0
             4         3          0 ];

    %��������նȾ���gK,����նȾ���gK
    node_number=length(gNode(:,1));
    element_number=length(gElement(:,1));
    gK=zeros(node_number*3,node_number*3);
    gM=zeros(node_number*3,node_number*3);
    gVout=zeros(node_number*3,pztn);
    [M,Kn,Ku]=controlrod();
    % ���㵥Ԫ�նȾ���k������װ�����������
    for ie=1:1:element_number
       [Kuu,ke,me]=ElementMatrixpzt(ie,m,gnum,mn,mtop,M,Kn,Ku);
       AssembleMatrix(ie,ke,me);
       if ie==1||ie==3||ie==5||ie==7
          Vout(ie,Kuu);
          cont=cont+1;
       end
    end
% ����Լ������
gK=BoundaryConditions(gBC,gK);  %�����������նȾ���
gM=BoundaryConditions(gBC,gM);  %�������������������
gVout=BC2(gBC,gVout);  %�������������������
% �����Ƶ�ʼ�����
[v,d]=eig(gK,gM);  
tempd=diag(d);
[nd,sortindex]=sort(tempd);
v=v(:,sortindex);
save para.mat gK gM gVout v;
save gnode.mat gNode gBC;
mode_number=1:20;    %��ע��ģ̬��
frequency(mode_number)=sqrt(nd(mode_number))/2/pi;
frequency = frequency';
toc;
    