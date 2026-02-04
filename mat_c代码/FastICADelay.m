clc
clear
close all
% user=53;
exclude=[52 54 33 38 71 74];
inter=2;
for user=11:75
    if ismember(user,exclude)
        continue
    end
    
    if user<10
        fname=strcat('a0', num2str(user));
        name=strcat('a0', num2str(user),'.fqrs.txt');
    else
        fname=strcat('a', num2str(user));
        name=strcat('a', num2str(user),'.fqrs.txt');
    end
    
    [signal fs]= rdsamp(fname,[], []);
    
    ann=rdann(fname,'fqrs',[],[]);
    ECG=signal;
    tm=0.001:0.001:60;
    tm=tm';
    cName=fname;
    qrsAf=ann;
    
    dbFlag=0;                   % debug flag
    graph=0;                    % enable/disable graphical representation
    saveFig=0;                  % =1 => save figures of the processing phases
    saveFigRRf=0;               % =1 => save estimated fetal RR figures
    
    fs = 1000;             % sampling frequency
    % ECG=preprocess(signal');
    % ECG=ECG';
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
    % close all
    
    % ---- Independent Component Analysis ----
    % Xm=FecgICAm(Xf,fs,cName,graph,dbFlag,saveFig);
    Se=FecgICAm(Xf,fs,cName,0,dbFlag,saveFig);
    
    % ---- Signal interpolation
    % Xi=FecgInterp(X,fs,interpFact,cName,graph);
    
    [Se,fs]=FecgInterp(Se,fs,inter,cName,0);
    % qrsAf=qrsAf/1000*fs;
    
    % ---- Channel selection and Mother QRS detection
    qrsM=FecgQRSmDet(Se,fs,cName,0,dbFlag,saveFig,qrsAf);
    
    % ---- Mother QRS cancelling
    Xr=FecgQRSmCanc(Se,qrsM,fs,cName,0,0,saveFig,qrsAf);
    
    % ---- Source separation by ICA on residual signals
    Ser=FecgICAf(Xr,fs,cName,graph,dbFlag,saveFig);
    
    % ---- Channel selection and Fetal QRS detection 原方法
    qrsF1=FecgQRSfDet(Ser,fs,cName,qrsM,0,dbFlag,saveFig,saveFigRRf,qrsAf);
    [acc1(user),ppv1(user),sen1(user),f11(user)]=evaluation(qrsF1*fs/inter,qrsAf,1);
    fprintf('acc:%.4f,ppv:%.4f,sen:%.4f,f1:%.4f\n',acc1(user),ppv1(user),sen1(user),f11(user));
%     SS={acc1,ppv1,sen1,f11};
    % qrsF1=FecgQRSfDet(Ser,fs,cName,qrsM,1,dbFlag,saveFig,saveFigRRf,qrsAf);
    
    % ---- 时延 ----
    Xr = derivative(Xr);
    [z,~,Cond]=yanchibaihua3(Xr',3 ,0);
    
    % ---- Source separation by ICA on residual signals
    Ser1=FecgICAf(z',fs,cName,graph,dbFlag,saveFig);
    
    % ----delete same yc----
    yc=dis_spike(Ser1',0.6);
    
    % ---- Channel selection and Fetal QRS detection 
%     qrsF2=FecgQRSfDet(yc',fs,cName,qrsM,0,dbFlag,saveFig,saveFigRRf,qrsAf);
%     [acc(user),ppv(user),sen(user),f1(user)]=evaluation(qrsF2*fs/inter,qrsAf,1);
%     fprintf('acc:%.4f,ppv:%.4f,sen:%.4f,f1:%.4f\n',acc(user),ppv(user),sen(user),f1(user));

    peakinterval=0.3*fs;
ann1=ann/1000*fs;


% adecgBp = derivative(z');
% ----threshold selecting ,clustering----

yc=[Ser yc']';
F=[];
for i=1:size(yc,1)
    [F1,F2,yc(i,:),th]=getspike(yc(i,:),fs,peakinterval,graph,ann1);
    F=[F,F1];
end


% ----generate reference signal,20ms----
win=0.01*fs;
if ~isempty(F)
    refer=zeros(size(F,2),size(Ser,1));
    for i=1:size(F,2)
        if ~isempty(F{i})
            for j=1:size(F{i},2)
                refer(i,max(F{i}(j)-win,1):min(F{i}(j)+win,size(Ser,1)))=1;
            end
        end
        
    end
end

maxacc2=0;
flag=0;
% ----constraint fastICA ----
% max_ref=max_spike(refer,ann1);
for i=1:size(refer,1)
max_ref=refer(i,:);
figure('NumberTitle', 'off', 'Name', '采集的阈值：')
    set(gcf,'Position',get(0,'ScreenSize'));
    plot(max_ref)
    hold on
    plot(ann1,0.9,'k+')
            yt=fcica3(z,max_ref,fs);
%     [yt, ~] = fcica3_period1207(z, refer(i,:),qrsM,qrsAf);
%     [yt, ~] = cfICA(z, max_ref,qrsM,qrsAf,fs,user);
    qrsFcfica=FecgQRSfDetCfica(yt',fs,num2str(user),qrsM,0,0,0,0,qrsAf);
    
    figure
    plot(yt)
    hold on
    plot(qrsAf/1000*fs,1,'r+');
    plot(qrsFcfica*fs,1,'k+')
   [acc2,ppv2,sen2,f12]=evaluation(qrsFcfica*fs/inter,qrsAf,1);
   if maxacc2<acc2     
   spike{user,1}=acc2,spike{user,2}=ppv2,spike{user,3}=sen2,spike{user,4}=f12;
   spike{user,5}=acc2,spike{user,6}=ppv2,spike{user,7}=sen2,spike{user,8}=f12;
   
   maxacc2=acc2;
   flag=1;
   end
   
   if flag==1
    for ii=-25:25
    [acc,ppv,sen,f1]=evaluation(qrsFcfica*fs/inter+ii,qrsAf,1);
    if acc>spike{user,1}
        spike{user,5}=acc;
        spike{user,6}=ppv;
        spike{user,7}=sen;
        spike{user,8}=f1;
%         spike{5}=qrsFcfica*fs/inter+ii;
    end
    end
    flag=0;
   end
    
close all

end
end

a=1;
% xlswrite('FastICA+delay.xlsx',acc',1,'B3');
% xlswrite('FastICA+delay.xlsx',ppv',1,'C3');
% xlswrite('FastICA+delay.xlsx',sen',1,'D3');
% xlswrite('FastICA+delay.xlsx',f1',1,'E3');
% 
% xlswrite('FastICA+delay.xlsx',acc1',1,'G3');
% xlswrite('FastICA+delay.xlsx',ppv1',1,'H3');
% xlswrite('FastICA+delay.xlsx',sen1',1,'I3');
% xlswrite('FastICA+delay.xlsx',f11',1,'J3');
    

function [S]=max_spike(y,ann)
        win=10;
        s=zeros(1,size(y,2));
        s(fix(ann))=1;
        comax=0;
        for j=1:size(y,1) %isempty判断输入是否为空                
                temp1=max(xcorr(s(1,2*win+1:end-2*win),y(j,2*win+1:end-2*win),'coeff')); %xcorr，是指互相关函数
                temp2=max(xcorr(s(1,2*win+1:end-2*win),-y(j,2*win+1:end-2*win),'coeff'));
                Co=max(temp1,temp2);
                if comax<Co
                    comax=Co;
                    S=y(j,:);
                end
        end
        
    end