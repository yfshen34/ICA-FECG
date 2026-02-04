clc
clear
close all

user=40;
inter=2;
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

% [Se,fs]=FecgInterp(Se,fs,inter,cName,0);
% qrsAf=qrsAf/1000*fs;

% ---- Channel selection and Mother QRS detection
qrsM=FecgQRSmDet(Se,fs,cName,0,dbFlag,saveFig,qrsAf);

% ---- Mother QRS cancelling
Xr=FecgQRSmCanc(Se,qrsM,fs,cName,0,0,saveFig,qrsAf);
load('a40_result.mat');
[XrF,XcF]=FecgQRSfCanc(Xr,fix(qrsff),fs,cName,0,0,saveFig,[]);

Xc=Se-Xr;
left=10*fs;
last=15;
aecg=Se(left:left+fs*last,2);
mecg=Xc(left:left+fs*last,2);
res=Xr(left:left+fs*last,2);
fecg=XcF(left:left+fs*last,2);

figure
set(gcf,'Position',get(0,'ScreenSize'));

subplot(3,1,1), plot(aecg,'linewidth',1.5);
        wgmi1= min(aecg) ;
        wgma1= max(aecg) ;
        ylim([wgmi1, wgma1]);
        set(gca,'Visible','off');
        
        
        subplot(3,1,2), plot(res,'linewidth',1.5);
        hold on
        plot(mecg,':r','linewidth',1.5);
        hold off
        ylim([wgmi1, wgma1]);
        set(gca,'Visible','off');       
        
        subplot(3,1,3), plot(fecg,'linewidth',1.5);
        ylim([wgmi1, wgma1]);
        hold on
        qrs=fix(qrsff(qrsff>=left ));
        qrs=qrs(qrs<=left+fs*last);
        qrs=qrs-left;
       
        [~,s]=findpeaks(fecg,'MINPEAKDISTANCE',350);
        plot(s,fecg(s),'or');
        hold off
       set(gca,'Visible','off');
        set(gcf,'Color','white');
set(gca,'Visible','off');
set(gcf,'Position', [100 100 2040 440]);%左，下，高度，宽度

       saveFig=1;
progname=mfilename;
if(saveFig), figFmt='png';
        figPath=fullfile('../Figure/',progname);
        if(~exist(figPath,'dir')), mkdir(figPath); end
        figName=fullfile(figPath,[cName,'svd']);
        print(gcf, ['-d',figFmt],figName);
        %    saveas(gcf, figName,'fig');
end
       