function AssembleMatrix(ie,ke,me)
%  �ѵ�Ԫ���󼯳ɵ��������
%  �������:
%      ie  --- ��Ԫ��
%      ke  --- ��Ԫ�նȾ���
%      me  --- ��Ԫ��������
%  ����ֵ:
%      ��
    global gElement gK gM
    for i=1:1:2   %��ֵȡ��Ԫ�����Ľڵ���
        for j=1:1:2    %��ֵȡ��Ԫ�����Ľڵ���
            for p=1:1:3   %��ֵȡ�ڵ�����ɶ���
                for q=1:1:3   %��ֵȡ�ڵ�����ɶ���
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
