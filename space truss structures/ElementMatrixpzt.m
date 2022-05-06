function [Kuu,k,ma] = ElementMatrixpzt(ie,m,gnum,mn,mtop,M,Kn,Ku)
global gNode gElement gMaterial
     n=gElement(ie,3);
     E=gMaterial(n,1);
     A=gMaterial(n,2); 
     density=gMaterial(n,3); 
     xi=gNode(gElement(ie,1),1); 
     yi=gNode(gElement(ie,1),2); 
     zi=gNode(gElement(ie,1),3);
     xj=gNode(gElement(ie,2),1);
     yj=gNode(gElement(ie,2),2);
     zj=gNode(gElement(ie,2),3);
     L=sqrt((xj-xi)^2+(yj-yi)^2+(zj-zi)^2); 
     lx=(xj-xi)/L; 
     mx=(yj-yi)/L;
     nx=(zj-zi)/L; % 为避免旋转奇异性，通过符号运算，直接给出变换后的表达式
     T=zeros(6);
     T(1,1)=lx;T(1,2)=mx;T(1,3)=nx;
     T(4,4)=lx;T(4,5)=mx;T(4,6)=nx;
     Kuu=T'*Ku;
     k=E*A/L*[ lx*lx  lx*mx  lx*nx  -lx*lx  -lx*mx  -lx*nx
                 mx*lx  mx*mx  mx*nx  -mx*lx  -mx*mx  -mx*nx
                 nx*lx  nx*mx  nx*nx  -nx*lx  -nx*mx  -nx*nx
                 -lx*lx  -lx*mx  -lx*nx  lx*lx  lx*mx  lx*nx
                 -mx*lx  -mx*mx  -mx*nx  mx*lx  mx*mx  mx*nx
                 -nx*lx  -nx*mx  -nx*nx  nx*lx  nx*mx  nx*nx ];
     me=density*A*L/6*[2 0 0 1 0 0;0 2 0 0 1 0;0 0 2 0 0 1;1 0 0 2 0 0
                0 1 0 0 2 0;0 0 1 0 0 2 ];
     if ie==1||ie==3||ie==5||ie==7
     k=Kn*[ lx*lx  lx*mx  lx*nx  -lx*lx  -lx*mx  -lx*nx
                 mx*lx  mx*mx  mx*nx  -mx*lx  -mx*mx  -mx*nx
                 nx*lx  nx*mx  nx*nx  -nx*lx  -nx*mx  -nx*nx
                 -lx*lx  -lx*mx  -lx*nx  lx*lx  lx*mx  lx*nx
                 -mx*lx  -mx*mx  -mx*nx  mx*lx  mx*mx  mx*nx
                 -nx*lx  -nx*mx  -nx*nx  nx*lx  nx*mx  nx*nx ];  
     me=M;
     end
%      m = me/2*eye(6);   %集中质量矩阵
if n==1&&ie/gnum>(m-1)
     ma=me+[0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 mtop 0 0
                0 0 0 0 mtop 0;0 0 0 0 0 mtop];%顶端质量附件杆
elseif n==1
     ma=me+[0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 mn 0 0
                0 0 0 0 mn 0;0 0 0 0 0 mn];%节点质量附加杆
else
     ma=me;%无附加质量杆
end  
end