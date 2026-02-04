function [c,ceq] = nlcon(x)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
c=25-x(1)*x(2)*x(3)*x(4);
ceq=sum(x.^2)-40;
end

