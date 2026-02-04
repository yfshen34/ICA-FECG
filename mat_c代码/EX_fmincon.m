% min func
objective=@(x) x(1)*x(4)*(x(1)+x(2)+x(3))+x(3);
 x0=[1,5,5,1];

disp(['initial objective:' num2str(objective(x0))])
A=[];
b=[];
Aeq=[];
Beq=[];
lb=1.0*ones(4);
ub=5.0*ones(4);

nonlincon=@nlcon;

[x,fvl,ef,output,lambda]=fmincon(objective,x0,A,b,Aeq,Beq,lb,ub,nonlincon);
output.iterations
disp(x)
disp(['final objective:' num2str(objective(x))])
[c ceq]=nlcon(x)
