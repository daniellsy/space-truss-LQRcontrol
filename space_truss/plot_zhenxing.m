% 绘制振型图
m = length(v)/3;
j=1; % 第几阶
k = 2;  % 放大系数
n = length(gNode(:,1))/4; % 共4*n个节点

for i=1:m/4
    x(4*(i-1)+1:4*(i-1)+4) = dx*i;
end

for i=1:m
    y(i) = k*v(3*(i-1)+2,j);
    z(i) = k*v(3*(i-1)+3,j);
end
% x = v(1:64,1);
% y = v(64+1:64*2,1);
% z = v(64*2+1:64*3,1);
xx = gNode(:,1);
y1=[zeros(1,4),y];
yy = gNode(:,2)+y1';
z1=[zeros(1,4),z];
zz = gNode(:,3)+z1';

% 邻接矩阵
A1 = [0 1 0 1;
      1 0 1 0;
      0 1 0 1;
      1 0 1 0];
B1 = [1 0 0 1;
      1 1 0 0;
      0 1 1 0;
      0 0 1 1];
B1 = kron(eye((n-1)),B1);
B1 = [zeros(4*(n-1),4),B1];
B1 = [B1;zeros(4,4*n)];
A = kron(eye(n),A1)+B1+B1';
            
coord1 = zeros(length(xx),3);
for i=1:length(xx)
    coord1(i,1)=xx(i);
    coord1(i,2)=yy(i);
    coord1(i,3)=zz(i);
end
figure(1)
gplot3(A,coord1,'color','black')
axis equal
axis([0 6.5 -1.5 1.9 -1.5 1.9]);
str = sprintf('第%1d阶振型图',j);
title(str);grid on
xlabel('x轴');ylabel('y轴');zlabel('z轴')
% plot3(xx,yy,zz);