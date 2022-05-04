function AssembleMatrix(ie,ke,me)
%  把单元矩阵集成到整体矩阵
%  输入参数:
%      ie  --- 单元号
%      ke  --- 单元刚度矩阵
%      me  --- 单元质量矩阵
%  返回值:
%      无
    global gElement gK gM
    for i=1:1:2   %终值取单元包含的节点数
        for j=1:1:2    %终值取单元包含的节点数
            for p=1:1:3   %终值取节点的自由度数
                for q=1:1:3   %终值取节点的自由度数
                    m = (i-1)*3+p ;
                    n = (j-1)*3+q ;
                    M = (gElement(ie,i)-1)*3+p ;
                    N = (gElement(ie,j)-1)*3+q ;
                    gK(M,N) = gK(M,N) + ke(m,n) ;
                    gM(M,N) = gM(M,N) + me(m,n) ;
                end
            end
        end
    end
end
