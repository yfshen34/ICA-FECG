function [f,gradf] = bigtoleft(x)
% 目标函数bigtoleft 辅助函数。
%   此处显示详细说明
w2=x0{1};
xw=x0{2};
refer=x0{3};

gaosi=randn(1,size(xw,2));
Ggs=mean(GG(gaosi));
f=-(GG());


    function G=GG(x)
        G=log(cosh(x));
        % G=x.*exp(-x.^2/2);
    end

    function G=g(x)
        G=tanh(x);
        % G=x.*exp(-x.^2/2);
    end

    function G=gd(x)
        G=1-(tanh(x)).^2;
        %  G=(1-x.^2).*exp(-x.^2/2);
    end
end

