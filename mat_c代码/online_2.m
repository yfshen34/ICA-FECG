clc
clearvars -except Qm Qf wm wf

% clear all
close all

user=5;
inter=2;
if user<10
    fname=strcat('a0', num2str(user));
    name=strcat('a0', num2str(user),'.fqrs.txt');
else
    fname=strcat('a', num2str(user));
    name=strcat('a', num2str(user),'.fqrs.txt');
end

[signal,fs]= rdsamp(fname,[], []);

ann=rdann(fname,'fqrs',[],[]);
ECG1=signal;
qrsf = [];
step = 1000;
window = 10000;
Xm0=[];
% 
imax = floor((60000 - window)/step)
for i=0:imax
    ECG=ECG1(step*i+1:step*i+window,:);
    tm=0.001:0.001:60;
    tm=tm';
    cName=fname;
    qrsAf=ann;
    
    dbFlag=0;                   % debug flag
    graph=0;                    % enable/disable graphical representation
    saveFig=0;                  % =1 => save figures of the processing phases
    saveFigRRf=0;               % =1 => save estimated fetal RR figures
    
    fs = 1000;             % sampling frequency
    % ---- Artifact canceling ----
    % 去尖峰噪声
    % X=FecgFecgImpArtCanc(ECG,fs,cName,graph,dbFlag);
    X=FecgImpArtCanc(ECG,fs,cName,0,saveFig);
    % close all
    
    % ---- detrending  ----低通滤波
    % Xd=FecgDetrFilt(X,fs,cName,graph,dbFlag);
    Xd=FecgDetrFilt(X,fs,cName,0,saveFig);
    % close all
    
    % ---- Power line interference removal by notch filtering ----
    % Xf=FecgNotchFilt(Xd,fs,cName,graph,dbFlag);
    Xf=FecgNotchFilt(Xd,fs,cName,0,saveFig);
    
    Se=Qm*Xf';
    Se=wm*Se;
    
    
    % ---- Signal interpolation
    [Se,fs]=FecgInterp(Se',fs,inter,cName,0);
    
    % ---- Channel selection and Mother QRS detection
    qrsM=FecgQRSmDet(Se,fs,cName,0,dbFlag,saveFig,qrsAf);
    
    % ---- Mother QRS cancelling
  %  Xr=FecgQRSmCanc(Se,qrsM,fs,cName,1,1,saveFig,qrsAf);
     [Xr,Xm]=FecgQRSmCanc(X,qrsM,fs,cName,1,1,saveFig,qrsAf);

    Xm1=Xm(501:550,1);
    Xm0=[Xm0;Xm1];
    % Xm1=Xm(:,indexm);
     for i=1:size(Xm,2)
        a=diff(Xm(:,i));
        % much=(length(a(a<1000))>=sum(a(a<1000))/fs*3);%发放个数的限制
        a((a>500))=[];
        b=Xm(qrsM,i);
        coef_s1(i)=std(a)/mean(a);
        coef_w1(i)=std(b)/abs(mean(b));
    end
    [cs1,idsm]=sort(coef_s1);%%s interval
      [cw1,idw]=sort(coef_w1);%% w amp
    sw=coef_s1+0.3*coef_w1;
    [csw,idsw]=sort(sw);
    Xm1=Xm(:,idsw(1));
    figure;
    % plot(Xm(:,2));
    % stackedplot(Xm);


  Xr = derivative(Xr);
    
    Xr1=Xr';
    [Nsignal, Nsample] = size(Xr');
    R=1;
    xt = zeros((2*R+1)*Nsignal,Nsample);
    for k= 1:Nsignal
        %xt(R*(i-1)+1:R*i,:)=toeplitz([x(i,1);zeros(R-1,1)],x(i,:)); toeplitz(x)用向量x生成一个对称的托普利兹矩阵
        xt((2*R+1)*(k-1)+1:(2*R+1)*k,:)=toeplitz([Xr1(k,(R+1):-1:1),zeros(1,R)],[Xr1(k,(R+1):end),zeros(1,R)]);%xt为扩展矩阵
    end
    xt=xt-mean(xt,2)*ones(1,Nsample);%将xt中心化 ,求每行的均值
    
    
    y=Qf*xt;
    
    Ser1=wf*y;
    Ser1=Ser1';
    
    stackedplot(Ser1); 


    qrsFcfica=FecgQRSfDet(Ser1,fs,cName,qrsM,graph,dbFlag,0,saveFigRRf,qrsAf);
    qrsF=qrsFcfica*1000+i*step;

    if (i==0)
        qrsAf1 = qrsAf(find(qrsAf>=step*i+1 & qrsAf<=step*i+window));
        qrsF1  = qrsF(find(qrsF>=step*i+1 & qrsF<=step*i+window));
        % err1 = qrsF1 - qrsAf1;
        % err1_mean = mean(abs(err1));
        qrsf = [qrsf;qrsF1];
    else 
        qrsAf1 = qrsAf(find(qrsAf>=step*(i-1)+window+1 & qrsAf<=step*i+window));
        qrsF1  = qrsF(find(qrsF>=step*(i-1)+window+1 & qrsF<=step*i+window));
        % err1 = qrsF1 - qrsAf1;
        % err1_mean = mean(abs(err1));
        qrsf = [qrsf;qrsF1];
    end
    % err=ans-qrsAf(44:65,1);
    % mean(abs(err));
    close all;
end
%err = qrsf - qrsAf;
%err_mean = mean(abs(err));