function [k,m] = ElementMatrix(ie,mn)
global gNode gElement gMaterial
     n= gElement(ie, 3);
     E = gMaterial( n, 1 ) ;
     A = gMaterial( n, 2 ) ; 
     density = gMaterial( n, 3 ) ; 
     xi = gNode( gElement( ie, 1 ), 1 ) ; 
     yi = gNode( gElement( ie, 1 ), 2 ) ; 
     zi = gNode( gElement( ie, 1 ), 3 ) ;
     xj = gNode( gElement( ie, 2 ), 1 ) ;
     yj = gNode( gElement( ie, 2 ), 2 ) ;
     zj = gNode( gElement( ie, 2 ), 3 ) ;
     L = sqrt((xj-xi)^2 + (yj-yi)^2 + (zj-zi)^2 ); 
     lx = abs((xj-xi)/L) ; 
     mx = abs((yj-yi)/L) ;
     nx = abs((zj-zi)/L) ;    % 为避免旋转奇异性，通过符号运算，直接给出变换后的表达式
     k = E*A/L*[ lx*lx  lx*mx  lx*nx  -lx*lx  -lx*mx  -lx*nx
                 mx*lx  mx*mx  mx*nx  -mx*lx  -mx*mx  -mx*nx
                 nx*lx  nx*mx  nx*nx  -nx*lx  -nx*mx  -nx*nx
                 -lx*lx  -lx*mx  -lx*nx  lx*lx  lx*mx  lx*nx
                 -mx*lx  -mx*mx  -mx*nx  mx*lx  mx*mx  mx*nx
                 -nx*lx  -nx*mx  -nx*nx  nx*lx  nx*mx  nx*nx ];
     me = density*A*L;
%      m = me/2*eye(6);   %集中质量矩阵
if n==1
     m = me/6*[ 2 0 0 1 0 0;0 2 0 0 1 0;0 0 2 0 0 1;1 0 0 2 0 0
                0 1 0 0 2 0;0 0 1 0 0 2 ]+...
                [ mn 0 0 0 0 0;0 mn 0 0 0 0;0 0 mn 0 0 0;0 0 0 0 0 0
                0 0 0 0 0 0;0 0 0 0 0 0 ]; % 一致质量矩阵
else
    m = me/6*[ 2 0 0 1 0 0;0 2 0 0 1 0;0 0 2 0 0 1;1 0 0 2 0 0
                0 1 0 0 2 0;0 0 1 0 0 2 ];
end
     %质量矩阵不需要变换（变换后也不发生改变）      
end