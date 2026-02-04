function [y, w] = cfICA(X, refer,qrsM,qrsAf,fs,user)
if nargin<5
    fs=1000;
end
inter=fs/1000;
spike={0,0,0,0,[]};

ref=refer/sqrt(mean(refer.^2));
fprintf('Starting cICA for extracting the desired source signal ..\n');
[ICnum, IClen]=size(X);

w0=mean(X(:,ref>0),2);
w0=w0/norm(w0);
w = w0;
oldw = w;

gamma1=1e1;
gamma2=1e1;

mu1 = 0;
mu2 = 0;


OverValue=0.0001;
maxIter = 200;

flag = 1;
loop = 1;
decrease=1;

% output at current iteration
    y = w' * X;


% compute the autocorrelation matrix Rxx
thr1 = 1 ;
r0 = xcorr(y);
    r0(1,floor(length(r0)/2)+1)=0;
    r0=r0(1,floor(length(r0)/2)+1:end);
    thr2=mean(power(r0,2));
    % 周期性的缩放因子
    scalar=thr2;
    RRchangeold=1;
    
% figure
while (flag == 1)
    % Update γ.
    gamma1=gamma1*1.02;
    gamma2=gamma2*1.02;
    
  
tic
    sumgd2=zeros(ICnum,1);
    r=xcorr(y);
    r(1,floor(length(r)/2)+1)=0;
    r=r(1,floor(length(r)/2)+1:end);
    
    if mu2==0
        gdg2=0;
    else
    for k=1:IClen-1
        if(k==1) 
            disp('进入大循环'); 
        end
        %     (Tk + TkT)*y 
        tkt = circshift(y',-k);
        tkt(end-k+1:end)=0;
        tk = circshift(y',k);
        tk(1:k)=0;
        T=tk+tkt;
        sumgd2=sumgd2+X*T*r(k+1);
    end
    gdg2=-2*sumgd2/IClen;
    end
   toc
    disp(['gdg21运行时间: ',num2str(toc)]);
    w =  ( X * tanh(y)')/IClen + mu1 *( X * ref')/IClen - mu2 *gdg2/scalar;
    w = w/norm(w);
    
    % output at current iteration
    y = w' * X;
    
    
    %     thr = threshold;
    g1 = thr1 - mean(y.*ref);
    if g1>0
        thr1 = thr1 *0.95;
    end
    
    g2=(thr2 - mean(power(r,2)))/scalar;
    if g2<0 
        thr2 = thr2 *1.03;
    end
    
    % update of the parameter mu1 
    mu1 = max(0, mu1 + gamma1 * g1);
    mu2 = max(0, mu2 + gamma2 * (g2-0.03));
      
    wchange = abs(abs(w'*oldw)-1);
    % fprintf('No.%d iteration: acc:%.4f,ppv:%.4f,sen:%.4f,f1:%.4f\n',loop, spike{1},spike{2},spike{3},spike{4});
    % fprintf('No.%d iteration: change in w is %g\n',loop, wchange);
    if wchange < OverValue
        fprintf('Converged after %d iteration\n',loop);
        flag = 0;
    end
    
    if loop >= maxIter
        fprintf('After %d iteration, still cannot convergent.\n',loop);
        flag = 0;
    end
    if decrease==1 && loop >=10
        maxIter=maxIter*10;
        decrease=0;
    end
    oldw = w;
    loop = loop + 1;
    close all
end

% output
y = w'* X;

% qrsFcfica=FecgQRSfDetCfica(y',fs,'02',qrsM,1,0,0,0,qrsAf);
%     figure
%     plot(y)
%     hold on
%     plot(qrsAf/1000*fs,1,'r+');
%     plot(qrsFcfica*fs,1,'k+')
% 
%     maxacc=0;
%         maxppv=0;
%         maxsen=0;
%         maxf1=0;
% for ii=-25:25
%     [acc,ppv,sen,f1]=evaluation(qrsFcfica*fs/inter+ii,qrsAf,1);
%     if acc>maxacc
%         maxacc=acc;
%         maxppv=ppv;
%         maxsen=sen;
%         maxf1=f1;
%     end
% end
% fprintf('acc:%.4f,ppv:%.4f,sen:%.4f,f1:%.4f\n',maxacc,maxppv,maxsen,maxf1);

%%%%%%%%%%%%%%%%%%%%%%%%心率计算%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% counter=1;
% left=1;
% last=5;
% qrs1=qrsFcfica;
% qrsFcfica1=qrs1;
% 
% qrs1=qrsFcfica1(qrs1>=left/fs);
% qrsFcfica1=qrsFcfica1(qrsFcfica1>=left/fs);
% qrs1=qrsFcfica1(qrs1<=(left+fs*last)/fs+last/fs);
% qrs1=qrs1-left/fs;
% qrs_i_raw=qrsFcfica1(1:length(qrs1)+counter );
% 
% temp = zeros(1,counter);
% % qrs_i_rand_start = 1;
% for jj=1:length(qrs1)
% for i = 1:counter
%     temp(i) = qrs_i_raw(jj+i)- qrs_i_raw(jj+i-1);
% %     temp(i) = t(qrs_i_raw(i+1))- t(qrs_i_raw(i));
% end
% % RR_Tavg is the average RR interval over a "counter" amount of samples
% % The unit is "second". It should be below 1 second. Normally around 0.5s
% RR_Tavg = sum(temp)/counter;
% 
% % Fetal Heart Rate per minute, physiological FHR should be around 120 to
% % 150 beats per minute
% FHR_per_min(jj,1) = 60/RR_Tavg;
% end
% 
% % mean(FHR_per_min)
% % std(FHR_per_min)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% fprintf('End of cICA algorithm !\n');



end
