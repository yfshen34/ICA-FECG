clc
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
step = 10;
window = 1000;
% 
% ECG = load("E:\project1\algo/algo_test_200hz_5000.txt");
%ECG = load("E:\project1\ECG_mat2c-dev-integration_android\data\filter_200Hz_35000.txt");
% ECG = load("E:\project1\ECG_mat2c-dev-integration_android\data\algo_test_200hz_5000.txt");
ECG=ECG1;
ECG = ECG(1001:13000,:);
fs = 200;
imax = floor((4000 - window)/step);
for i=0:0
% for i=8
    % ECG1=ECG(step*i+1:step*i+window,:);
    tm=0.001:0.001:60;
    tm=tm';
    cName=fname;
    qrsAf=ann;
    
    dbFlag=0;                   % debug flag
    graph=0;                    % enable/disable graphical representation
    saveFig=0;                  % =1 => save figures of the processing phases
    saveFigRRf=0;               % =1 => save estimated fetal RR figures
    
    fs = 200;             % sampling frequency
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
    
    while(max(Xf,[],"all")<=50 || max(Xf,[],"all")>=200)
        if max(Xf)<=50
            Xf = Xf*2;
        elseif max(Xf)>=200
            Xf = Xf/2;
        end
    end

    Se=Qm*Xf';
    Se=wm*Se;
    
    
    % ---- Signal interpolation
    [Se1,fs]=FecgInterp(Se',fs,inter,cName,0);
    
    % ---- Channel selection and Mother QRS detection
    qrsM=FecgQRSmDet(Se1,fs,cName,0,dbFlag,saveFig,qrsAf);
    
    % ---- Mother QRS cancelling
    [Xr2,Xe]=FecgQRSmCanc(Se1,qrsM,fs,cName,1,1,saveFig,qrsAf);
    %%其中Xe是母心

    Xr = derivative(Xr2,fs);
    
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

    qrsFcfica=FecgQRSfDet(Ser1,fs,cName,qrsM,graph,dbFlag,0,saveFigRRf,qrsAf);
    % qrsF=qrsFcfica*fs/2+i*step;
    qrsF=qrsFcfica*fs;
    [Xres,resi] = FecgQRSfCanc(Ser1,qrsF,fs,cName,graph,dbFlag,saveFig);
    

    % if (i==0)
    %     qrsAf1 = qrsAf(find(qrsAf>=step*i+1 & qrsAf<=step*i+window));
    %     qrsF1  = qrsF(find(qrsF>=step*i+1 & qrsF<=step*i+window));
        % err1 = qrsF1 - qrsAf1;
        % err1_mean = mean(abs(err1));
    %     qrsf = [qrsf;qrsF1];
    % else 
    %     qrsAf1 = qrsAf(find(qrsAf>=step*(i-1)+window+1 & qrsAf<=step*i+window));
    %     qrsF1  = qrsF(find(qrsF>=step*(i-1)+window+1 & qrsF<=step*i+window));
        % err1 = qrsF1 - qrsAf1;
        % err1_mean = mean(abs(err1));
        % qrsf = [qrsf;qrsF1];
    % end
    % err=ans-qrsAf(44:65,1);
    % mean(abs(err));
    close all;
    fprintf('in 5 s,qrsM=%5d\n',qrsM);
    fprintf('in 5 s,qrsM.size()=%4d\n', size(qrsM,1));
    fprintf('in 5 s,qrsF=%5d\n',qrsF);
    fprintf('in 5 s,qrsF.size()=%4d\n', size(qrsF,1));
end
% qrsAf_window = qrsAf(find(qrsAf<=step*i+window));
% err = qrsf - qrsAf_window;
% err_mean = mean(abs(err));