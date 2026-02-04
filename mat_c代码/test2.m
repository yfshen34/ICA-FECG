clc
clear 
close all

user=7;
inter=2;
if user<10
    fname=strcat('a0', num2str(user)); %字符连接
    name=strcat('a0', num2str(user),'.fqrs.txt');
else
    fname=strcat('a', num2str(user));
    name=strcat('a', num2str(user),'.fqrs.txt');
end

[signal fs]= rdsamp(fname,[], []); %读信号

ann=rdann(fname,'fqrs',[],[]); %读标注
ECG=signal;
% ECG=ECG(1:5000,:);
% ECG=emg(7:10,100000:150000)';
% ECG=resample(ECG,1,5);
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

% ---- Power line interference removal by notch filtering ----%陷波滤波
% Xf=FecgNotchFilt(Xd,fs,cName,graph,dbFlag);
Xf=FecgNotchFilt(Xd,fs,cName,0,saveFig);
% close all
% Xm=FecgICAm(Xf,fs,cName,graph,dbFlag,saveFig); %ICA独立成分分析
[Se,Qm,wm]=FecgICAm(Xf,fs,cName,0,dbFlag,saveFig);

% ---- Signal interpolation %信号插值
[Se,fs]=FecgInterp(Se,fs,inter,cName,0);

% ---- Channel selection and Mother QRS detection %检测母体心电
qrsM=FecgQRSmDet(Se,fs,cName,0,dbFlag,saveFig,qrsAf);

% ---- Mother QRS cancelling %SVD消除母体心电
Xr1=FecgQRSmCanc(Se,qrsM,fs,cName,1,1,saveFig,qrsAf);

Xr = derivative(Xr1,fs); %差分滤波


 Fs=1000;
adaptiveFF.profile = 'cooling';
adaptiveFF.tau_const = Inf; % unit: samples
% pars for cooling ff
adaptiveFF.gamma = 0.6;
adaptiveFF.lambda_0 = 0.995; 
% pars for adaptive ff
adaptiveFF.decayRateAlpha = 0.02;
adaptiveFF.upperBoundBeta = 1e-3;
adaptiveFF.transBandWidthGamma = 1;
adaptiveFF.transBandCenter = 5;
adaptiveFF.lambdaInitial = 0.1;
 if strcmp(adaptiveFF.profile,'cooling') || strcmp(adaptiveFF.profile,'constant')
        adaptiveFF.lambda_const  = 1-exp(-1/(adaptiveFF.tau_const)); % steady state constant lambda
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% emgdata=syfu10;
evalConvergence.profile = true;
evalConvergence.leakyAvgDelta = 0.01; % Leaky average value (delta) for computing non-stationarity index (NSI). NSI = norm(Rn), where Rn = (1-delta)*Rn + delta*(I-yf^T).
evalConvergence.leakyAvgDeltaVar = 1e-3; % Leaky average value (delta) for computing variance of source activation. Var = (1-delta)*Var + delta*variance.
onlineWhitening = false;
X_tmp=Xf;
nChs=size(X_tmp,2);
nPts=8;
n=size(X_tmp,1);
 state.lambda_k      = zeros(1,40);   % readout lambda
    state.minNonStatIdx = []; 
    state.counter       = 0; % time index counter, used to keep track of time for computing lambda
    state.Rn = [];
    state.nonStatIdx = [];
    state.icaweights = eye(nChs);
    state.icasphere = eye(nChs);
    state.icasphere = 2.0*inv(sqrtm(double(cov(X_tmp))));
     state.kurtsign      = ones(nChs,1) > 0;
        %%%%%%%%%%%%%%%%%%%%%%%
data =  state.icasphere*X_tmp';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% i=1;
% res=[];
% x1=data;
% while (i-1)*8<=size(x1,2)
%     dataRange=(i-1)*40+1:mni(i*40,n);
%     xx = x1(:,dataRange);
%    state = dynamicOrica(xx, state, dataRange, adaptiveFF, evalConvergence);
%    onys=state.icaweights*xx;
% %     if i == 1
% %         res = [res,onys];
% %     else
% %         res = [res,onys(:,end-399:end)];
% %     end
%    res = [res,onys];
%     i = i+1;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step=40;
% numBlock=floor(n/step);
% res=[];
% cres=[];
% miu=0;
% sigma=1e1;ksi=1;
% cw=zeros(size(40,1),1);
% for i = 0 : numBlock-1
%        
%         dataRange = 1 + floor(i*step) : min(n, floor((i+1)*step));
%         xx = data(:,dataRange);
% %         [state,miu,sigma,ksi,cw] = ocica2(xx, state, dataRange, adaptiveFF, evalConvergence,miu,sigma,ksi,cw);
%           state = dynamicOrica(xx, state, dataRange, adaptiveFF, evalConvergence);
% 
%         onys=state.icaweights'*xx;
% %         conys=cw'*xx;
%         res=[res,onys];
% %         cres = [cres,conys];
% 
% end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


step=200;
numBlock=floor(n/step);
res=[];
cres=[];
icaEMG1res=[];
huanyuanres=[];
  sigma=1e1;ksi=1; miu=0;
 
% miu=0;
% sigma=1e1;ksi=1;
% % cw=zeros(size(40,1),1);
lamdba=1;
xw=yanchibaihua3(data,1,0.0);

%% 

for i = 0 : numBlock-1
       
        dataRange = 1 + floor(i*step) : min(n, floor((i+1)*step));
        xx = data(:,dataRange);
%         [state, miu,sigma,ksi,cw] = ocica2(xx, state, dataRange, adaptiveFF, evalConvergence, miu,sigma,ksi,cw);
       state = dynamicOrica(xx, state, dataRange, adaptiveFF, evalConvergence);
% xishu=state.icaweights*state.icasphere;
        onys=state.icaweights*data(:,dataRange);
        
       %增加cICA约束 
            res = [res,onys];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    zuizhongjieguo=state.icaweights*data;
%% 


 PlotEMG(res, Fs, [], [], 'Independent Components');
%   PlotEMG(huanyuanres, Fs, [], [], 'Independent Components');
%   plot(cres)
%   baoluo=envelope(cres,50,'rms');
% t=(1:length(res))/Fs;
% figure;plot(t,cres);box off;ax1=gca;ax1.YAxis.Visible = 'off';
% xlim([0 63])
xlabel('Time(s)')
 set (0,'defaultfigurecolor','w') 