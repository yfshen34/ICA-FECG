
% close all

%  user=5;
% inter=2;
% % for user = 1:75
% %     try
% 
% if user<10
%     fname=strcat('a0', num2str(user));
%     name=strcat('a0', num2str(user),'.fqrs.txt');
% else
%     fname=strcat('a', num2str(user));
%     name=strcat('a', num2str(user),'.fqrs.txt');
% end
% 
% [signal fs]= rdsamp(fname,[], []);
% 
% ann=rdann(fname,'fqrs',[],[]);
% load('b1515.mat');
% signal=data15(20001:56000,:);

% filename = '1119_2.csv';       % 你的 CSV 文件
% T = readtable(filename);       % 读取为 table（带表头）
% 
% signal = T{:,:};
% 
% save('1119-2.mat', 'signal');


ECG=signal(10001:49000,:);
% ECG(:,1)=ECG(:,2)+ECG(:,3);
% ECG= ECG(1:20000,:);
 fs=200;

tm=0.001:0.001:60;
tm=tm'; %----画波形的？
fname=0;
cName=fname;
ann=0;
qrsAf=ann;
cName=0;
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
Xd=FecgDetrFilt(X,fs,cName,graph,dbFlag);
% Xd=FecgDetrFilt(X,fs,cName,0,saveFig);
% close all

% ---- Power line interference removal by notch filtering ----
Xf=FecgNotchFilt(Xd,fs,cName,graph,dbFlag);
% Xf=FecgNotchFilt(Xd,fs,cName,0,saveFig);
% close all

% ---- Independent Component Analysis ----
% Xm=FecgICAm(Xf,fs,cName,graph,dbFlag,saveFig);
%Se 源信号
%Qm 白化矩阵，wm 解混矩阵，这两个参数离线保存，在线用到
[Se,Qm,wm]=FecgICAm(Xf,fs,cName,0,dbFlag,saveFig);
% [Se,Qm,wm]=FecgICAm(Xf,fs,0,0,dbFlag,saveFig);


% ---- Signal interpolation
% [Se,fs]=FecgInterp(Se,fs,inter,cName,0);


% ---- Channel selection and Mother QRS detection
qrsM=FecgQRSmDet(Se,fs,cName,0,dbFlag,saveFig,qrsAf);


% ---- Mother QRS cancelling
%Xr1=FecgQRSmCanc(Se,qrsM,fs,cName,1,1,saveFig,qrsAf);
[Xr1,Xm]=FecgQRSmCanc(Se,qrsM,fs,0,1,1,0,qrsAf);
% [Xr1,Xm]=FecgQRSmCanc(Xf,qrsM,fs,0,1,1,0,qrsAf);

% Xr = derivative(Xr1,fs);
[z,Qf,~]=yanchibaihua3(Xr1',1 ,0);

% ---- Source separation by ICA on residual signals
% Ser=FecgICAf(Xr,fs,cName,graph,dbFlag,saveFig);
[Ser1,~,wf]=FecgICAf(z',fs,cName,graph,dbFlag,saveFig);

qrsFcfica=FecgQRSfDet(Ser1,fs,cName,qrsM,graph,dbFlag,0,saveFigRRf,qrsAf);
qrsf=floor(qrsFcfica*fs);
qrsMica=qrsM/fs;
Xf1=FecgQRSfCanc(Xr1,qrsf,fs,cName,1,1,saveFig,qrsAf);
RRm=diff(qrsMica);
RRf=diff(qrsFcfica);
mhr=60/mean(RRm);
fhr=60/mean(RRf);

valid_idx_m = (RRm >= 0.3) & (RRm <= 2);  % 过滤异常值（生理合理范围：0.3~2秒）
RRm_clean = RRm(valid_idx_m);             % 有效RR间期

if numel(RRm_clean) >= 2
    mean_RRm = mean(RRm_clean);   % 有效RR间期的均值
    std_RRm = std(RRm_clean);     % 有效RR间期的标准差
    cv_RRm = std_RRm / mean_RRm;  % 母QRS的CV
else
    cv_RRm = NaN;                 % 有效数据不足时设为NaN
    warning('母QRS的RR间期有效数据不足（<2个），无法计算CV');
end

valid_idx_f = (RRf >= 0.3) & (RRf <= 2);  % 过滤异常值
RRf_clean = RRf(valid_idx_f);             % 有效RR间期

if numel(RRf_clean) >= 2
    mean_RRf = mean(RRf_clean);   % 有效RR间期的均值
    std_RRf = std(RRf_clean);     % 有效RR间期的标准差
    cv_RRf = std_RRf / mean_RRf;  % 胎儿QRS的CV
else
    cv_RRf = NaN;                 % 有效数据不足时设为NaN
    warning('胎儿QRS的RR间期有效数据不足（<2个），无法计算CV');
end

% --- 母心 MHR（母亲）---
rr_min_m = 0.3;   % 最低 RR 间期（秒） ≈ 200 bpm
rr_max_m = 2.0;   % 最高 RR 间期（秒） ≈ 30 bpm
valid_idx_m = (RRm >= rr_min_m) & (RRm <= rr_max_m);

% --- 胎儿 FHR（胎儿）---
rr_min_f = 0.3;   % 最低 RR 间期（秒） ≈ 200 bpm
rr_max_f = 1.0;   % 最高 RR 间期（秒） ≈ 60 bpm （推荐上限，可设为 1.2s）
% rr_max_f = 1.2; % 可选：如果你想更宽松一点，设为 1.2s（≈50 bpm）

valid_idx_f = (RRf >= rr_min_f) & (RRf <= rr_max_f);

RRm_clean = RRm(valid_idx_m);
if numel(RRm_clean) >= 2
    mean_RRm = mean(RRm_clean);
    std_RRm = std(RRm_clean);
    cv_RRm = std_RRm / mean_RRm;
else
    cv_RRm = NaN;
    warning('母心 RR 间期有效数据不足（<2个），无法计算 CV');
end

% 胎儿：筛选后计算
RRf_clean = RRf(valid_idx_f);
if numel(RRf_clean) >= 2
    mean_RRf = mean(RRf_clean);
    std_RRf = std(RRf_clean);
    cv_RRf = std_RRf / mean_RRf;
else
    cv_RRf = NaN;
    warning('胎儿 RR 间期有效数据不足（<2个），无法计算 CV');
end

% --- 1. 离群 RR 间期（更严格范围，精细筛选）
rr_min_m_tight = 0.6;
rr_max_m_tight = 1.2; %母心范围50-100
is_outlier_m_tight = (RRm_clean < rr_min_m_tight) | (RRm_clean > rr_max_m_tight);
outlier_ratio_m_tight = sum(is_outlier_m_tight) / length(RRm_clean);

rr_min_f_tight = 0.375;
rr_max_f_tight = 0.6; %胎心范围100-160
is_outlier_f_tight = (RRf_clean < rr_min_f_tight) | (RRf_clean > rr_max_f_tight);
outlier_ratio_f_tight = sum(is_outlier_f_tight) / length(RRf_clean);

fprintf('离群母心RR（母心范围50-100）比例: %.2f%%\n', outlier_ratio_m_tight * 100);
fprintf('离群胎儿RR（胎心范围100-160）比例: %.2f%%\n', outlier_ratio_f_tight * 100);

% --- 2. 相邻 RR 间期变化幅度过大（动态阈值 30%）
mean_rr_m = mean(RRm_clean);
mean_rr_f = mean(RRf_clean);
threshold_change = 0.40;

delta_rr_m = abs(diff(RRm_clean));
large_change_m = delta_rr_m > (threshold_change * mean_rr_m);
large_change_ratio_m = sum(large_change_m) / length(delta_rr_m);

delta_rr_f = abs(diff(RRf_clean));
large_change_f = delta_rr_f > (threshold_change * mean_rr_f);
large_change_ratio_f = sum(large_change_f) / length(delta_rr_f);

fprintf('【RR变化】母心相邻RR间期变化过大比例(>%.0f%% 均值): %.2f%%\n', threshold_change*100, large_change_ratio_m * 100);
fprintf('【RR变化】胎儿相邻RR间期变化过大比例(>%.0f%% 均值): %.2f%%\n', threshold_change*100, large_change_ratio_f * 100);
% ica_filename = ['ica_params_', fname, '.mat'];
% save(ica_filename, 'Qm', 'wm', 'Qf', 'wf','user','fname');
% catch ME
% warning('Error processing %s: %s', fname, ME.message);
%         continue;
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%心率计算%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
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
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%心率计算%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% counter=1;
% left=1;
% last=10;
% qrs1=qrsM/2000;
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
% for jk=1:length(qrs1)
% for i1 = 1:counter
%     temp(i1) = qrs_i_raw(jk+i1)- qrs_i_raw(jk+i1-1);
% %     temp(i) = t(qrs_i_raw(i+1))- t(qrs_i_raw(i));
% end
% % RR_Tavg is the average RR interval over a "counter" amount of samples
% % The unit is "second". It should be below 1 second. Normally around 0.5s
% RR_Tavg = sum(temp)/counter;
% 
% % Fetal Heart Rate per minute, physiological FHR should be around 120 to
% % 150 beats per minute
% FHR_per_minM(jk,1) = 60/RR_Tavg;
% end



% % ----delete same yc----
% yc=dis_spike(Ser1',0.6);
% 
% 
% peakinterval=0.3*fs;
% debug=0;
% ann1=ann/1000*fs;
% 
% 
% % ----threshold selecting ,clustering----
% 
% yc=[Ser yc']';
% % yc=Ser';
% F=[];
% for i=1:size(yc,1)
%     [F1,F2,yc(i,:),th]=getspike(yc(i,:),fs,peakinterval,debug,ann1);
%     F=[F,F1];
% end
% 
% 
% 
% % ----generate reference signal,20ms----
% win=0.01*fs;
% if ~isempty(F)
%     refer=zeros(size(F,2),size(Ser,1));
%     for i=:1size(F,2)
%         if ~isempty(F{i})
%             for j=1:size(F{i},2)
%                 refer(i,max(F{i}(j)-win,1):min(F{i}(j)+win,size(Ser,1)))=1;
%             end
%         end
% 
%     end
% end
% 
% % ----constraint fastICA ----
% i=1;
% while i<=size(refer,1)
%     %         yc(i,:)=fcica3(xw,refer(i,:));
%     figure('NumberTitle', 'off', 'Name', '采集的阈值：')
%     set(gcf,'Position',get(0,'ScreenSize'));
%     plot(refer(i,:))
%     hold on
%     plot(ann1,0.9,'k+')
%     % de=input('Please input the threshold(or press "enter" to auto setting or input "0" to delete it):','s');
%     % de = str2num(de);
%     % if de~=1
%     %     fprintf('跳过第%d个Spike\n', i);
%     %     i=i+1;
%     %     continue
%     % end
%     [yt, ~] = cfICA(z, refer(i,:),qrsM,qrsAf,fs,user);
% 
%     qrsFcfica=FecgQRSfDet(yt',fs,cName,qrsM,graph,dbFlag,0,saveFigRRf,qrsAf);
% 
% 
%     ycmax = max(yt);
%     figure;set(gcf,'Position',get(0,'ScreenSize'));
%     hold on;
%     plot(yt);
% 
%     title([num2str(i),'/',num2str(size(yt,1))]);
%     plot(ann1,1,'r+');
%     plot(qrsFcfica*fs,1,'k+');
%     evaluation(qrsFcfica*fs/inter,qrsAf,1)
%     figResize(0, 1, 1, .35);
%     i=i+1;
% 
% 
% 
% close all
% end


