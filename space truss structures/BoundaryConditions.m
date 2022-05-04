function matrix= BoundaryConditions(gBC,Matrix)   
%  处理位移约束条件（划行划列法）
%  输入参数：
%     gBC -- 节点位移约束矩阵
%     Matrix -- 要缩减的矩阵
%  返回值：
%     matrix  -- 缩减之后的矩阵

bc_number = length(gBC(:,1));
list = zeros(bc_number,1);    %存储要删去的行列号
for i = 1:bc_number
    list(i) = (gBC(i,1)-1)*3+gBC(i,2); % *节点自由度数
end
    Matrix(list,:) = [];
    Matrix(:,list) = [];
    matrix = Matrix;
end