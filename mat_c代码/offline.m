clc
% clear 
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

[signal fs]= rdsamp(fname,[], []);

ann=rdann(fname,'fqrs',[],[]);
ECG=signal;
% ECG=resample(ECG,1,5);
% ECG=load("G:\data\filter_200Hz_35000.txt");
ECG=ECG(1:10000,:);
% ECG=emg(7:10,100000:150000)';

fs=1000;
% ECG=ECG(1:10000,:);
tm=0.001:0.001:60;
tm=tm';
cName=fname;
qrsAf=ann;

dbFlag=0;                   % debug flag
graph=0;                    % enable/disable graphical representation
saveFig=0;                  % =1 => save figures of the processing phases
saveFigRRf=0;               % =1 => save estimated fetal RR figures

% fs = 1000;             % sampling frequency
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
[Se,Qm,wm]=FecgICAm(Xf,fs,cName,0,dbFlag,saveFig);

% ---- Signal interpolation
[Se,fs]=FecgInterp(Se,fs,inter,cName,0);

% ---- Channel selection and Mother QRS detection
qrsM=FecgQRSmDet(Se,fs,cName,0,dbFlag,saveFig,qrsAf);

% ---- Mother QRS cancelling
Xr1=FecgQRSmCanc(Se,qrsM,fs,cName,1,1,saveFig,qrsAf);

Xr = derivative(Xr1,fs);
[z,Qf,~]=yanchibaihua3(Xr',1 ,0);
z1=z;
% ---- Source separation by ICA on residual signals
% Ser=FecgICAf(Xr,fs,cName,graph,dbFlag,saveFig);

[Ser,~,wf]=FecgICAf(z',fs,cName,graph,dbFlag,saveFig);


peakinterval=0.3*fs;
debug=0;
ann1=ann/1000*fs;


% F=[];
% yc=Ser';
% for i=1:size(yc,1)
%     [F1,F2,yc(i,:),th]=getspike(yc(i,:),fs,peakinterval,debug,ann1);
%     F=[F,F1];
% end



% ----delete same yc----
yc=dis_spike(Ser',0.6);



% ----threshold selecting ,clustering----

% yc=Ser';
F=[];
for i=1:size(yc,1)
    [F1,F2,yc(i,:),th]=getspike(yc(i,:),fs,peakinterval,debug,ann1);
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

% ----constraint fastICA ----
i=1;
qrsFcfica=cell(1,size(refer,1));
wt=zeros(size(refer,1),size(z,1));
while i<=size(refer,1)

    [yt, wtt] = cfICA(z, refer(i,:),qrsM,qrsAf,fs,user);
    wt(i,:)=wtt;
    % qrsFcfica{i}= FecgQRSfDet(yt',fs,cName,qrsM,graph,dbFlag,0,saveFigRRf,qrsAf);

    i=i+1;

close all
end
yt1=wt*z;
 qrsf=FecgQRSfDet(yt1',fs,cName,qrsM,graph,dbFlag,0,saveFigRRf,qrsAf);
 qrsf=floor(qrsf*fs);
%  [Xr11,resi]=FecgQRSfCanc(Ser1,qrsf,fs,cName,1,1,saveFig,qrsAf);

