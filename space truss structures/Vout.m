function Vout(ie,Kuu)
%  �ѵ�Ԫ���󼯳ɵ��������
%  �������:
%      ie  --- ��Ԫ��
%      ke  --- ��Ԫ�նȾ���
%      me  --- ��Ԫ��������
%  ����ֵ:
%      ��
    global gElement gVout cont
    for i=1:1:2   %��ֵȡ��Ԫ�����Ľڵ���
            for p=1:1:3   %��ֵȡ�ڵ�����ɶ���
                    m=(i-1)*3+p;
                    M=(gElement(ie,i)-1)*3+p;
                    gVout(M,cont)=gVout(M,cont)+Kuu(m);
            end
    end
end


