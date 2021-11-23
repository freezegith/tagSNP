%清空环境变量  
clear;
clc; 
%导入SNP数据
SNP = importdata('17964357-163.haps');
%取出tagSNP训练数据

a1=84;a2=85;a3=10;a4=30;a5=16;
ysize=112;

tagSNP_train = [SNP(1:1600,a1),SNP(1:1600,a2),SNP(1:1600,a3),SNP(1:1600,a4),SNP(1:1600,a5)];
%取出tagSNP测试数据
tagSNP_test = [SNP(1601:2000,a1),SNP(1601:2000,a2),SNP(1601:2000,a3),SNP(1601:2000,a4),SNP(1601:2000,a5)];
%加入一列元素全为“1”的列向量 x1，构造 X
x1((1:1600), 1) = 1;
X = [x1, tagSNP_train];
%预测准确率
accuracy(1:ysize, 1) = 0;
B(1:6, 1:ysize) = 0;
Y_result(1:400, 1:ysize) = 0;
for x = 1:ysize
    if x==a1 || x==a2 || x==a3 || x==a4 || x==a5 
        continue;
    end
    %取出 Y的训练数据
    Y = [SNP(1:1600, x)];
    %取X、Y矩阵的行数、列数
    [X_r, X_c] = size(X);
    [Y_r, Y_c] = size(Y);
    %b为参数，bint为回归系数的区间估计，r为残差，rint为置信区间，stats用于回归模型检验 
    [b,bint,r,rint,stats]=regress(Y,X); 
    B(:,x) = b;
    [b_r, b_c] = size(b);
    %加入一列元素全为“1”的列向量 x2，构造 X_test
    x2((1:400), 1) = 1;
    X_test = [x2, tagSNP_test];
    %取出 Y的测试数据
    Y_test = [SNP(1601:2000, x)];
    %计算预测结果 y
    y((1:400), 1) = 0;
    %预测正确的个数
    num_true = 0;
    %预测总个数
    allNum = 400;

    for i = 1:400
        for j = 1:b_r
            y(i) = y(i) + b(j) * X_test(i, j);
        end
        Y_result(i,x) = y(i);
        %取近似值
        if y(i) < 0.3
            y(i) = 0;
        elseif y(i) > 0.7
            y(i) = 1;
        else
            allNum = allNum - 1;
            continue;
        end
        %验证正确性
        if y(i) == Y_test(i)
            num_true = num_true + 1;
        end
    end
    %disp('预测准确率:');
    accuracy(x) = num_true / allNum;
end
xlswrite('accuracy.xls', accuracy);