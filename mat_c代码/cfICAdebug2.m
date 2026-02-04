% function [y, w] = cfICA(X, refer)

clc
clear
close all

load('a02bestyc.mat')
yc=yc';
load('2residuebaihua.mat')

S=searchspike_in(yc,1.7,300,1000,1,ann);
evaluation(S,ann);
% load('2residue.mat')
% load('2fcica.mat')
a=find(S==1);
for i=1:length(a)
    S(a(i)-20:a(i)+20)=1;
end
% S(S==0)=-0.1;    

figure
plot(S)
hold on
plot(ann,1,'r+')
% plot(yc)

%  load('2residue.mat')
% X=z;
% X=double(xw);
% X=whiten(xw);
x=[residue(1,:);residue(2,:); residue(3,:)];

qrs=find(Spike==1);
for i=1:length(qrs)
    x(:,qrs(i)-15:qrs(i)+15)=0;
end

[z,~,Cond]=yanchibaihua3(x,1 ,0);
X=z;

ref=S/sqrt(mean(S.^2));
fprintf('Starting cICA for extracting the desired source signal ..\n');
[ICnum, IClen]=size(X);

minth1=mean(S.*ref)

w0=mean(X(:,ref>0),2);
w0=w0/norm(w0);
w = w0;
oldw = w;

gamma1=1e1;
gamma2=1e1;

mu1 = 0;
mu2 = 0;

OverValue=1e-5;
maxIter = 200;

flag = 1;
loop = 1;

% output at current iteration
    y = w' * X;

% compute the autocorrelation matrix Rxx
% Rxx = X(:,[1:end]) * X(:,[1:end])' / IClen;
thr1 = 0.99 ;
r = xcorr(y);
    r(1,floor(length(r)/2)+1)=0;
    r=r(1,floor(length(r)/2)+1:end);
    thr2=mean(power(r,2));
    % 周期性的缩放因子
    scalar=thr2;

    
figure
 plot(y)
 
 [yc, ~] = fcica3_period1207(X,ref);
 plot(yc)
hold on
plot(ann,1,'r+')
hold off
S=searchspike_in(yc,-1,320,1000,1,ann);
evaluation(S,ann);
while (flag == 1)
    % Update γ.
    gamma1=gamma1*1.02;
    gamma2=gamma2*1.02;
    
     % calculate the first order deviation of the Lagarange function
%     std_y = std(y);                          % standard deviation
%     v_Gaus = normrnd(0, std_y, 1, IClen);    % Gaussian signal with the same mean and variance
%     rou = mean( log(cosh(y)) - log(cosh(v_Gaus)) );
    %     L1 = sign(rou) * ( X * tanh(y)')/IClen - mu1 * ( X * (y - ref)')/IClen ...
    %         - lambda * ( X * y')/IClen;
    
    
    % related to the second order deviation of the Lagarange function
    %     Sw = sign(rou) * mean(1-tanh(y).^2) - mu1 - lambda;
    
    % update of the weight vector
    %     w = w - learningRate * inv(Rxx) * L1 / Sw;
    %     w = w - learningRate * L1 / Sw;
    %     w =  ( X * tanh(y)')/IClen - mu1 *( X * ref')/IClen;

    %     rt=zeros(IClen,IClen);
%     r=zeros(1,IClen);   
tic
%     gdg2x=zeros(ICnum,IClen);
    sumgd2=zeros(ICnum,1);
    r=xcorr(y);
    r(1,floor(length(r)/2)+1)=0;
    r=r(1,floor(length(r)/2)+1:end);
    
    if mu2==0
        gdg2=0;
    else
%     gdg2x(:,1)=X*(y'+y');
    for k=1:IClen-1
        %     (Tk + TkT)*y 
        tkt = circshift(y',-k);
        tkt(end-k+1:end)=0;
        tk = circshift(y',k);
        tk(1:k)=0;
        T=tk+tkt;
%         gdg2x(:,k+1)=X*T;
        sumgd2=sumgd2+X*T*r(k+1);
        
                
%         r(k)=y*yk;
%         gdg2x(:,k)=X*yk;
    end
    gdg2=-2*sumgd2/IClen;
    
%     gdg21=-2*gdg2x*r'/IClen;
    end
   toc
    disp(['gdg21运行时间: ',num2str(toc)]);
%     gdg2=-2*X*rt*r'/IClen;
    w =  ( X * tanh(y)')/IClen + mu1 *( X * ref')/IClen - mu2 *gdg2/scalar;
    w = w/norm(w);
    
    % output at current iteration
    y = w' * X;
    
    
    %     thr = threshold;
    g1 = thr1 - mean(y.*ref);
    g2=(thr2 - mean(power(r,2)))/scalar;
    if g1>0 && g1>minth1
        thr1 = thr1 *0.95;
    elseif g2<0  && g1>minth1
        thr2 = thr2 *1.05;
        
    end
    
%     g2=(thr2 - mean(power(r,2)))/scalar;
%     if g2<0 
%         thr2 = thr2 *1.05;
%     end
    
    % update of the parameter mu1 
    %  g = mean( (y-ref).^2 ) - thr;    % corresponds to the inequality constraint
    mu1 = max(0, mu1 + gamma1 * g1);
    mu2 = max(0, mu2 + gamma2 * (g2-0.03));
    
    % update of the parameter lambda
    %     h = mean(y.^2) - 1;                    % corresponds to the equality constraint
    %     lambda = lambda + gamma * h;
    
    plot(y)
    pause(3)
    % decide whether the algorithm has converged or not
    wchange = abs(abs(w'*oldw)-1);
    fprintf('No.%d iteration: change in w is %g\n',loop, wchange);
    if wchange < OverValue
        fprintf('Converged after %d iteration\n',loop);
        flag = 0;
    end
    
    if loop >= maxIter
        fprintf('After %d iteration, still cannot convergent.\n',loop);
        flag = 0;
    end
    
    oldw = w;
    loop = loop + 1;
    
end

% output
y = w'* X;
g1
g2
% figure
% plot(y)
fprintf('End of cICA algorithm !\n');
S=searchspike_in(y,-1.1,350,1000,1,ann,4.4);
evaluation(S,ann);
figure
plot(y)
hold on
plot(ann,1,'r+')
hold off


% y1=abs(y);
y1=y;
qrs=find(Spike==1);
for i=1:length(qrs)
    y1(qrs(i)-20:qrs(i)+20)=0;
end
plot(y1)
hold on
plot(ann,1,'r+')
hold off
S=searchspike_in(y,-1,320,1000,1,ann,4.4);
evaluation(S,ann);