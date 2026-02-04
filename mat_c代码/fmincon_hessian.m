load('fcica.mat')
S=refer(3,:);
[channel,t]=size(xw);
%拉格朗日函数的 Hessian 矩阵
% 该函数首先计算 ?2f(x)。然后，它计算两个约束：Hessian ?2c1(x) 和 ?2c2(x)，将它们乘以对应的拉格朗日乘数 lambda.ineqnonlin(1) 和 lambda.ineqnonlin(2)，并将它们相加。
%bigtoleft 是一个目标函数
options = optimoptions('fmincon','Algorithm','interior-point',...
    "SpecifyConstraintGradient",true,"SpecifyObjectiveGradient",true,...
    'HessianFcn',@hessinterior);
w2=mean(xw(:,S>0),2);w2=w2/norm(w2);
x0{1} = w2;
x0{2}=xw;
x0{3}=refer;
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];

[x,fval,eflag,output] = fmincon(@bigtoleft,x0,...
           A,b,Aeq,beq,lb,ub,@twocone,options);
       
disp([output.funcCount,output.iterations])       
       
       