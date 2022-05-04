%    ����ƽ���ϵ������Ԫģ�����ݣ�
%        gNode ------- ��㶨��
%        gElement ---- ��Ԫ����
%        gMaterial --- ���϶��壬��������ģ�������Ľ���� 
%        gBC --------- Լ������ 
%        gK  --------- ����նȾ���
%        gM ---------  ������������
clear;clc;tic;
global gNode gElement gMaterial gBC  gK gM 
    % ����ڵ�͵�Ԫ
    %                         x����   y����
    gNode =  importdata('Input node coordinates.txt');    %  ��ȡ���ڵ�����
    %                         �ڵ��1   �ڵ��2   ���ϱ��
    gElement =  importdata('Input node and material numbers.txt'); %  ��ȡ����Ԫ�ڵ��
     
     % ��������
     %           ����ģ��   �����   �ܶ�  ���ɱ�
     gMaterial = [7.0e10   2.0e-4   2700   0.3       % ���� 1 (����
                  2.1e11   3.0e-4   7800   0.3];     % ���� 2���֣�
     
     % Լ������
     %     ����   ���ɶȺ�    Լ��ֵ 
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

    %��������նȾ���gK,����նȾ���gK
    node_number = length(gNode(:,1)) ;
    element_number = length(gElement(:,1)) ;
    gK = zeros( node_number * 3, node_number * 3 ) ;
    gM = zeros( node_number * 3, node_number * 3 ) ;
    
    % ���㵥Ԫ�նȾ���k������װ�����������
    for ie=1:1:element_number
       [ke,me] = ElementMatrix1(ie);
       AssembleMatrix(ie,ke,me);
    end
    
    % ����Լ������
    Ka=BoundaryConditions(gBC,gK);  %�����������նȾ���
    Ma=BoundaryConditions(gBC,gM);  %�������������������
    
% �����Ƶ�ʼ�����
[v,d]=eig(Ka,Ma);  
tempd=diag(d);
[nd,sortindex]=sort(tempd);
v=v(:,sortindex);
mode_number=1:20;    %��ע��ģ̬��
frequency(mode_number)=sqrt(nd(mode_number))/2/pi;
frequency = frequency';
toc;
    