function matrix= BoundaryConditions(gBC,Matrix)   
%  ����λ��Լ�����������л��з���
%  ���������
%     gBC -- �ڵ�λ��Լ������
%     Matrix -- Ҫ�����ľ���
%  ����ֵ��
%     matrix  -- ����֮��ľ���

bc_number = length(gBC(:,1));
list = zeros(bc_number,1);    %�洢Ҫɾȥ�����к�
for i = 1:bc_number
    list(i) = (gBC(i,1)-1)*3+gBC(i,2); % *�ڵ����ɶ���
end
    Matrix(list,:) = [];
    Matrix(:,list) = [];
    matrix = Matrix;
end