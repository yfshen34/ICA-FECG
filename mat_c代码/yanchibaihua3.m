function [y,Q,Cond] = yanchibaihua3(x, R,po)
% function xt= yanchibaihua(x, R)
%白化去相关
% 这里白化的意思是：
% 对原始信号X进行预处理 预处理包括去均值和白化(Whitening)，是通过对观测数据向量进行线性变换，使其均值为零，方差为1，去除各观测之间的相关性。
%
% 白化过程 简单而言就是 将信号或者噪声的协方差矩阵的对角化处理，不同的信号或信号形式其协方差矩阵一般而言是不一样的，那么我们从矩阵论的知识可以得到，两者要求的正交矩阵也是不一样的。
%R是the delay factor of constrained FastICA, the default value is 20 samples
%x为原始信号 Nsignal通道数
[Nsignal, Nsample] = size(x);



xt = zeros((2*R+1)*Nsignal,Nsample);
for k= 1:Nsignal
    %xt(R*(i-1)+1:R*i,:)=toeplitz([x(i,1);zeros(R-1,1)],x(i,:)); toeplitz(x)用向量x生成一个对称的托普利兹矩阵
    xt((2*R+1)*(k-1)+1:(2*R+1)*k,:)=toeplitz([x(k,(R+1):-1:1),zeros(1,R)],[x(k,(R+1):end),zeros(1,R)]);%xt为扩展矩阵
end
xt=xt-mean(xt,2)*ones(1,Nsample);%将xt中心化 ,求每行的均值
%  mxt = xt*xt'/Nsample;

% xt=x;

% x=x-mean(x,2)*ones(1,Nsample);

x_cov=cov(xt');                    % cov为求协方差的函数
%当Cond很大时，可认为x_cov为奇异矩阵，即行列式为0.
Cond=cond(x_cov);%矩阵A的条件数等于A的范数与A的逆的范数的乘积，即cond(A)=‖A‖·‖A^(-1)‖，对应矩阵的3种范数，相应地可以定义3种条件数。 函数 cond(A,1)、cond(A)或cond(A inf) 是判断矩阵病态与否的一种度量，条件数越大矩阵越病态。
[E,D]=eig(x_cov);  % 对信号矩阵的协方差函数进行特征值分解 E是特征向量 D是特征值(只有对角线有值)
% [E,D]=eig(mxt);
lamda=diag(D); %diag是(提取对角元素)
%将特征根及其特征向量按降序的方式重新排列
num=length(lamda);trun=fix(num*(1-po));
[w,j]=sort(lamda,'descend');E=E(:,j(1:trun)');D=diag(w(1:trun));%D排序好了，E每行顺序调换后和D顺序一致


Q=sqrt(D)\(E)'; % Q为白化矩阵
y=Q*xt;   %y为正交矩阵    MixedS_white为白化后的信号矩阵

% IsI=cov(y');                   %近似单位阵
end

