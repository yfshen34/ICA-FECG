function [y,xg,y0,miu,sigma,ksi]=fcica3_period1207(xw,S,qrsM,qrsAf)
% load('fcica.mat')
% S=yc(3,:);
%约束:q=(ksi^4-E{(yr)^4})/4
fs=1000;
% Set ? to a small value to control the learning precision.
% 迭代终点
thr=1e-5;

%是否结束迭代
flag=1;

% 迭代次数
c=40;

T=size(xw,2);
% reference signal
r=S/sqrt(mean(S.^2));

figure
while flag
    % 惩罚因子γ
    sigma=0.1;
    %时域约束阈值
    ksi=1;
    % Initialize wi with a random vector; set the norm to 1.
    %w1=rand(size(xw,1),1);w1=w1/norm(w1);
    w1=zeros(size(xw,1),1);
    w2=mean(xw(:,S>0),2);w2=w2/norm(w2);
    y0=w2'*xw;
    qrsFcfica=FecgQRSfDet(y0',1000,'02',qrsM,1,0,0,0,qrsAf);
    evaluation(qrsFcfica*fs,qrsAf,1)
    
    
    % For the i-th constrained IC: Initialize μ.
    miu=1;%the Lagrangian multiplier μ
    k=1;%迭代次数
    aa=abs(abs(w1'*w2)-1);
    % Repeat the following steps until the norm of Δwi is less than ? or
    % |1?||wi,kT *wi,k+1|||≤? (k is the iteration index).

    while aa>thr
        % Update γ.
        sigma=sigma*1.02;
        
        % Update wi using Eq. (15).
        w1=w2;
        y=w1'*xw;
        w2=1/T*xw*g(y)'+miu/T*xw*r';
        
        % Normalize wi using Eqs. (10) and (11).
        w2=w2/norm(w2);        
        
        
        %  Δwi is less than ksi， Adjust closeness threshold to gradually release the constraints.
        xg=mean((w2'*xw.*r));
        if xg<=ksi
            ksi=ksi*0.95;
        end
        
        % Update μ using Eq. (12a).        
        miu=max([0,miu+sigma*(ksi-xg)]);
        
        
        k=k+1;
        % 仍不收敛，放大收敛界限
        if k>c
            thr=thr*10;
            break;
        end
        aa=abs(abs(w1'*w2)-1);
    end
    %收敛成功，退出迭代
    flag=0;
    %需要在有限次数中收敛成功，若超过40次需要重新来过？
    if k>c
        c=c+10;
        flag=1;
    end
    %显示结果
    win=5;
    yy=w2'*xw;
    qrsFcfica=FecgQRSfDet(yy',1000,'02',qrsM,1,0,0,0,qrsAf);
    evaluation(qrsFcfica*fs,qrsAf,1)

end
w=w2;
y=w'*xw;
% figure;subplot(212);plot(y);subplot(211);plot(r);
%         function G=G(x)
%             %G=log(cosh(x));
%             G=x.*exp(-x.^2/2);
%         end

function G=g(x)
% G=tanh(x);
G=x.*exp(-x.^2/2);
end
%         function G=gd(x)
%             %G=1-(tanh(x)).^2;
%             G=(1-x.^2).*exp(-x.^2/2);
%         end
    end

