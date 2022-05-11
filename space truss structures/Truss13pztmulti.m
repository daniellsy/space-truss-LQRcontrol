%    定义平面杆系的有限元模型数据：
%        gNode ------- 结点定义
%        gElement ---- 单元定义
%        gMaterial --- 材料定义，包括弹性模量，梁的截面积 
%        gBC --------- 约束条件 
%        gK  --------- 整体刚度矩阵
%        
clear;clc;tic;
global gNode gElement gMaterial gpMaterial gBC gK gM gVout cont
cont=1;pztn=4;%压电主动杆单元数
% 材料性质
%             弹性模量   截面积   密度  泊松比
% gMaterial = [7.0e10   2.0e-4   2700   0.3       % 材料 1 (铝）
%              7.0e10   2.0e-4   2700   0.3       % 材料 2 (铝）
%              2.1e11   3.0e-4   7800   0.3];     % 材料 3（钢）
gMaterial = [7.27e10   6.283e-5   3100   0.3       % 横杆材料 1 (铝）
             7.27e10   6.283e-5   3100   0.3       % 竖杆材料 2 (铝）
             7.27e10   6.283e-5   3100   0.3];     % 斜杆材料 3（钢）
gpMaterial=[7600 7e-4 6.3e10 0.16 150 6.45e-10 1.5e-8 8000 6.283e-5 21e10 0.1];
%           roup    Ap     Ep    lp   np   d33     miup   rou    A       E    l
dx=0.3;%伸长方向上的一节的长度
dy=0.4;%宽度方向上一节的长度
m=21;%节数
mn=0;%接头质量
mtop=0;
gNode=zeros((m+1)*4,3);    
for i=1:m+1
      gNode((i-1)*4+1,1:3)=[dx*(i-1),0,0];  %节点坐标（x,y,z）
      gNode((i-1)*4+2,1:3)=[dx*(i-1),dy,0];  %节点坐标（x,y,z）
      gNode((i-1)*4+3,1:3)=[dx*(i-1),dy,dy];  %节点坐标（x,y,z）
      gNode((i-1)*4+4,1:3)=[dx*(i-1),0,dy];  %节点坐标（x,y,z）
end
 % 定义单元
gnum=13;%定义每个单元有几根杆
gElement=zeros(m*gnum,3);
for i=1:m
        gElement((i-1)*gnum+1,1:3)=[(i-1)*4+1,(i-1)*4+5,1];%横杆为1 竖杆为2 斜杆为3
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
     % 约束条件
     %     结点号   自由度号    约束值 
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

    %定义整体刚度矩阵gK,整体刚度矩阵gK
    node_number=length(gNode(:,1));
    element_number=length(gElement(:,1));
    gK=zeros(node_number*3,node_number*3);
    gM=zeros(node_number*3,node_number*3);
    gVout=zeros(node_number*3,pztn);
    [M,Kn,Ku]=controlrod();
    % 计算单元刚度矩阵k，并组装到整体矩阵中
    for ie=1:1:element_number
       [Kuu,ke,me]=ElementMatrixpzt(ie,m,gnum,mn,mtop,M,Kn,Ku);
       AssembleMatrix(ie,ke,me);
       if ie==1||ie==3||ie==5||ie==7
          Vout(ie,Kuu);
          cont=cont+1;
       end
    end
% 处理约束条件
gK=BoundaryConditions(gBC,gK);  %缩减后的整体刚度矩阵
gM=BoundaryConditions(gBC,gM);  %缩减后的整体质量矩阵
gVout=BC2(gBC,gVout);  %缩减后的整体质量矩阵
% 求固有频率及振型
[v,d]=eig(gK,gM);  
tempd=diag(d);
[nd,sortindex]=sort(tempd);
v=v(:,sortindex);
save para.mat gK gM gVout v;
save gnode.mat gNode gBC;
mode_number=1:20;    %关注的模态数
frequency(mode_number)=sqrt(nd(mode_number))/2/pi;
frequency = frequency';
toc;
    