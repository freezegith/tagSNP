%��ջ�������  
clear;
clc; 
%����SNP����
SNP = importdata('17964357-163.haps');
%ȡ��tagSNPѵ������

a1=84;a2=85;a3=10;a4=30;a5=16;
ysize=112;

tagSNP_train = [SNP(1:1600,a1),SNP(1:1600,a2),SNP(1:1600,a3),SNP(1:1600,a4),SNP(1:1600,a5)];
%ȡ��tagSNP��������
tagSNP_test = [SNP(1601:2000,a1),SNP(1601:2000,a2),SNP(1601:2000,a3),SNP(1601:2000,a4),SNP(1601:2000,a5)];
%����һ��Ԫ��ȫΪ��1���������� x1������ X
x1((1:1600), 1) = 1;
X = [x1, tagSNP_train];
%Ԥ��׼ȷ��
accuracy(1:ysize, 1) = 0;
B(1:6, 1:ysize) = 0;
Y_result(1:400, 1:ysize) = 0;
for x = 1:ysize
    if x==a1 || x==a2 || x==a3 || x==a4 || x==a5 
        continue;
    end
    %ȡ�� Y��ѵ������
    Y = [SNP(1:1600, x)];
    %ȡX��Y���������������
    [X_r, X_c] = size(X);
    [Y_r, Y_c] = size(Y);
    %bΪ������bintΪ�ع�ϵ����������ƣ�rΪ�вrintΪ�������䣬stats���ڻع�ģ�ͼ��� 
    [b,bint,r,rint,stats]=regress(Y,X); 
    B(:,x) = b;
    [b_r, b_c] = size(b);
    %����һ��Ԫ��ȫΪ��1���������� x2������ X_test
    x2((1:400), 1) = 1;
    X_test = [x2, tagSNP_test];
    %ȡ�� Y�Ĳ�������
    Y_test = [SNP(1601:2000, x)];
    %����Ԥ���� y
    y((1:400), 1) = 0;
    %Ԥ����ȷ�ĸ���
    num_true = 0;
    %Ԥ���ܸ���
    allNum = 400;

    for i = 1:400
        for j = 1:b_r
            y(i) = y(i) + b(j) * X_test(i, j);
        end
        Y_result(i,x) = y(i);
        %ȡ����ֵ
        if y(i) < 0.3
            y(i) = 0;
        elseif y(i) > 0.7
            y(i) = 1;
        else
            allNum = allNum - 1;
            continue;
        end
        %��֤��ȷ��
        if y(i) == Y_test(i)
            num_true = num_true + 1;
        end
    end
    %disp('Ԥ��׼ȷ��:');
    accuracy(x) = num_true / allNum;
end
xlswrite('accuracy.xls', accuracy);