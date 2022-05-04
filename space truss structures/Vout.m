function Vout(ie,Kuu)
%  把单元矩阵集成到整体矩阵
%  输入参数:
%      ie  --- 单元号
%      ke  --- 单元刚度矩阵
%      me  --- 单元质量矩阵
%  返回值:
%      无
    global gElement gVout cont
    for i=1:1:2   %终值取单元包含的节点数
            for p=1:1:3   %终值取节点的自由度数
                    m=(i-1)*3+p;
                    M=(gElement(ie,i)-1)*3+p;
                    gVout(M,cont)=gVout(M,cont)+Kuu(m);
            end
    end
end


