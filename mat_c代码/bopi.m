function [residual,huifu,huifu_ECG,W]=bopi(y,S,wavelength,p)%pwei1zehuatu
%residual:残差信号
%huifu：恢复得到的总信号
%huifu_ECG：恢复得到的ECG信号
%W：解混矩阵
a=randi(size(y,1),1,1);
if isempty(S)
    residual=y;
    W=[];
    huifu=y;
    return
end
[~,xt,~,huifu_ECG,C]=boxing(y(a,:),S,wavelength,p);
huifu=y*C'*xt;
residual=y-huifu;
W=C*y';

    function [W,xt,huifu,huifu_ECG,C]=boxing(x,S,wavelength,p)%x行向量。spike矩阵
        %% STA波形估计 LS
        % 动态确定波形长度
        [n,T]=size(S);
        yanchi=fix(0.45*wavelength);
        xt=zeros(n*wavelength,T);xt=single(xt);S=single(S);%xt是一个扩展矩阵
        
        for j=1:n
            %托普利兹矩阵
            %%%从yanchi+1往前取到1，因此toeplitz中第一个参数大小为wavelength，第二个参数大小即为S
            xt((j-1)*wavelength+1:j*wavelength,:)=toeplitz([S(j,(yanchi+1):-1:1), ...
                zeros(1,(wavelength-yanchi-1))], ...
                [S(j,(yanchi+1):end),zeros(1,yanchi)]);
        end
        C=(xt*xt')\xt;
        W=C*x';
        huifu=W'*xt;
        
        %             xls = pinv(x_in'*x_in)*x_in'*y_in;
        %             % Generate LAD coeff
        %             % Setup reformulated LP variables
        %             A = xt';
        %             b = x';
        %
        %             % LAD-TO-LP REFORMULATION with suppression trick
        %             len_x = size(A,2);
        %             c = [zeros(len_x,1);ones(size(b))];
        %             F = [A -eye(size(b,1)); -A -eye(size(b,1))];
        %             g = [b; -b];
        %             %F = [A -eye(size(b,1)); -A -eye(size(b,1));zeros(size(A)) -eye(size(b,1))];
        %             %g = [b ; -b; zeros(length(b),1)];
        %
        %             % Run the LP solver
        %             z = linprog(c,F,g);
        %             xlad = z(1,:);
        
        %在此处进行ECG绘制
        huifu_ECG=[];
        for j=1:n
            huifu_ECG=[huifu_ECG;W((j-1)*wavelength+1:j*wavelength,:)'*xt((j-1)*wavelength+1:j*wavelength,:)];
        end
        
        if (p==1)
            figure;set(gcf,'Position',get(0,'ScreenSize'));
            maxi=1.5*max(max(abs(x)),max(abs(huifu)));
            subplot(2,4,2:4);plot(x);hold on;
            for j=1:n
                plot(huifu_ECG(j,:)-j*maxi);hold on;
                %                     plot(huifu_ECG(:,:,j)-j*maxi);hold on;
            end
            plot(huifu-(j+1)*maxi);hold on;
            plot(x-huifu-(j+2)*maxi);hold on;
            subplot(2,4,6:8);
            for ii=1:size(S,1)
                plot(S(ii,:)-2*ii);hold on;
            end
            hold on;
            subplot(n,4,1:4:4*n-3);jg=max(abs(W));
            for j=1:n
                %                     plot(Ws(:,:,j)-1*jg*j,'linewidth',3);hold on;
                plot(W(1+wavelength*(j-1):wavelength*j)-2*jg*j,'linewidth',3);hold on;
            end
            pause(2);
        end
    end

end
