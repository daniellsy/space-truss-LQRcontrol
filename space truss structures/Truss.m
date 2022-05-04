%    定义平面杆系的有限元模型数据：
%        gNode ------- 结点定义
%        gElement ---- 单元定义
%        gMaterial --- 材料定义，包括弹性模量，梁的截面积 
%        gBC --------- 约束条件 
%        gK  --------- 整体刚度矩阵
%        gM ---------  整体质量矩阵
clear;clc;tic;
global gNode gElement gMaterial gBC  gK gM 
    % 定义节点和单元
    %                         x坐标   y坐标
    gNode =  importdata('Input node coordinates.txt');    %  读取各节点坐标
    %                         节点号1   节点号2   材料编号
    gElement =  importdata('Input node and material numbers.txt'); %  读取各单元节点号
     
     % 材料性质
     %           弹性模量   截面积   密度  泊松比
     gMaterial = [7.0e10   2.0e-4   2700   0.3       % 材料 1 (铝）
                  2.1e11   3.0e-4   7800   0.3];     % 材料 2（钢）
     
     % 约束条件
     %     结点号   自由度号    约束值 
     gBC = [ 1         1          0 
             1         2          0
             1         3          0
             5         1          0 
             5         2          0
             5         3          0
             9         1          0
             9         2          0
             9         3          0
             13        1          0
             13        2          0
             13        3          0 ];

    %定义整体刚度矩阵gK,整体刚度矩阵gK
    node_number = length(gNode(:,1)) ;
    element_number = length(gElement(:,1)) ;
    gK = zeros( node_number * 3, node_number * 3 ) ;
    gM = zeros( node_number * 3, node_number * 3 ) ;
    
    % 计算单元刚度矩阵k，并组装到整体矩阵中
    for ie=1:1:element_number
       [ke,me] = ElementMatrix1(ie);
       AssembleMatrix(ie,ke,me);
    end
    
    % 处理约束条件
    Ka=BoundaryConditions(gBC,gK);  %缩减后的整体刚度矩阵
    Ma=BoundaryConditions(gBC,gM);  %缩减后的整体质量矩阵
    
% 求固有频率及振型
[v,d]=eig(Ka,Ma);  
tempd=diag(d);
[nd,sortindex]=sort(tempd);
v=v(:,sortindex);
mode_number=1:20;    %关注的模态数
frequency(mode_number)=sqrt(nd(mode_number))/2/pi;
frequency = frequency';
toc;
    