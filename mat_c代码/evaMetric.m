clc
clear
close all

%合并
fs=1000;
ffid = fopen('RECORDS','r');
index=1;
while ~feof(ffid)
    tline{index} = fgetl(ffid);  %将每一行放到tline里
    index=index+1;
end
fclose(ffid);

i=1;

s=split(tline{i},'.');
% 导入文件
[hdr, record] = edfread(tline{i});

%导入注释
ann_name=strcat(s{1},'.txt');
ann =txtReader(ann_name);

for ii=1:size(ann,1)-1
    qrs(ii,1)=ann{ii+1,3};
end
qrs=qrs*fs;

cName=strcat('a',s{1},'_',int2str(i));
dbFlag=0;                   % debug flag
graph=1;                    % enable/disable graphical representation
saveFig=0;                  % =1 => save figures of the processing phases
saveFigRRf=0;               % =1 => save estimated fetal RR figures


qrsF=[];
for j=1:5
    qrs1=[];
    fname=strcat('a',s{1},'_',int2str(j),'_result.mat');
    
    load(fname);
    if j==1
        qrsF=qrsff(qrsff<=65*fs);
    elseif j==5
        qrs1=qrsff(qrsff<=60*fs);
        qrs1=qrs1(qrs1>=5*fs);
        qrs1=qrs1+(j-1)*60*fs;
        qrsF=[qrsF;qrs1];
    else
        qrs1=qrsff(qrsff<=65*fs);
        qrs1=qrs1(qrs1>=5*fs);
        qrs1=qrs1+(j-1)*60*fs;
        qrsF=[qrsF;qrs1];
    end
    
end
mother=[];
for k=1:2
    qrs1=[];
    fname=strcat('qrsM',s{1},'_',int2str(k),'.mat');
    load(fname);
    if k==1
        mother=qrsM(qrsM<=155*fs);
    else
        qrs1=qrsM(qrsM<=150*fs);
        qrs1=qrs1(qrs1>5*fs);
        qrs1=qrs1+(k-1)*150*fs;
        mother=[mother;qrs1];
    end
end
% qrsF=qrsF-15;
[acc,ppv,sen,f1]=evaluation(qrsF,qrs,1);


%-----------心率曲线---------------
counter=3;

%%%%%%%%%%%%%%%%%%%  MHR  %%%%%%%%%%%%%%%%%%%%%%%%%
qrsAf=mother/fs;
qrsFcfica=qrsAf;

qrs1=qrsFcfica(qrsFcfica<=105);
qrs1=qrsFcfica(qrs1>=55);
qrs1=qrs1-55;
qrs_i_raw=qrsFcfica(1:length(qrs1)+counter );

temp = zeros(1,counter);
% qrs_i_rand_start = 1;
for jj=1:length(qrs1)
for i = 1:counter
    temp(i) = qrs_i_raw(jj+i)- qrs_i_raw(jj+i-1);
%     temp(i) = t(qrs_i_raw(i+1))- t(qrs_i_raw(i));
end
% RR_Tavg is the average RR interval over a "counter" amount of samples
% The unit is "second". It should be below 1 second. Normally around 0.5s
RR_Tavg = sum(temp)/counter;

% Fetal Heart Rate per minute, physiological FHR should be around 120 to
% 150 beats per minute
FHR_per_min(jj,1) = 60/RR_Tavg;
end

mean(FHR_per_min)
std(FHR_per_min)

qrs1=[0;qrs1];
FHR_per_min=[FHR_per_min(1);FHR_per_min];
figure
subplot(4,1,1)
hold on
set(gcf,'color','w');
box off
% set(gca,'YColor','white')
ylabel('FHR(bpm)');
xlabel('Time(s)');
set(gca,'XTick',[0:10:50])

for i=2:length(qrs1)
    line([qrs1(i-1),qrs1(i-1)],[FHR_per_min(i-1),FHR_per_min(i)],'color','b','LineWidth',1.5);
    line([qrs1(i-1),qrs1(i)],[FHR_per_min(i),FHR_per_min(i)],'color','b','LineWidth',1.5);
end

ylim([70,100]);
set(gca,'YTick',[70:10:100])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qrsAf=qrs/fs;
qrsFcfica=qrsAf;
% qrs1=qrsFcfica(qrsFcfica<=60);
% 55001:105000,:
qrs1=qrsFcfica(qrsFcfica<=105);
qrs1=qrsFcfica(qrs1>=55);
qrs1=qrs1-55;
qrs_i_raw=qrsFcfica(1:length(qrs1)+counter );

temp = zeros(1,counter);
% qrs_i_rand_start = 1;
for jj=1:length(qrs1)
for i = 1:counter
    temp(i) = qrs_i_raw(jj+i)- qrs_i_raw(jj+i-1);
%     temp(i) = t(qrs_i_raw(i+1))- t(qrs_i_raw(i));
end
% RR_Tavg is the average RR interval over a "counter" amount of samples
% The unit is "second". It should be below 1 second. Normally around 0.5s
RR_Tavg = sum(temp)/counter;

% Fetal Heart Rate per minute, physiological FHR should be around 120 to
% 150 beats per minute
FHR_per_min(jj,1) = 60/RR_Tavg;
end

qrs1=[0;qrs1];
FHR_per_min=[FHR_per_min(1);FHR_per_min];
figure
subplot(4,1,1)
hold on
set(gcf,'color','w');
box off
% set(gca,'YColor','white')
ylabel('FHR(bpm)');
xlabel('Time(s)');
set(gca,'XTick',[0:10:50])
for i=2:length(qrs1)
    line([qrs1(i-1),qrs1(i-1)],[FHR_per_min(i-1),FHR_per_min(i)],'color','r','LineWidth',4);
    line([qrs1(i-1),qrs1(i)],[FHR_per_min(i),FHR_per_min(i)],'color','r','LineWidth',4);
end
ylim([120,150]);

qrsF=qrsF/fs;
qrsFcfica=qrsF;
% qrs1=qrsFcfica(qrsFcfica<=60);
% 55001:105000,:
qrs1=qrsFcfica(qrsFcfica<=105);
qrs1=qrsFcfica(qrs1>=55);
qrs1=qrs1-55;
qrs_i_raw=qrsFcfica(1:length(qrs1)+counter );

temp = zeros(1,counter);
% qrs_i_rand_start = 1;
for jj=1:length(qrs1)
for i = 1:counter
    temp(i) = qrs_i_raw(jj+i)- qrs_i_raw(jj+i-1);
%     temp(i) = t(qrs_i_raw(i+1))- t(qrs_i_raw(i));
end
% RR_Tavg is the average RR interval over a "counter" amount of samples
% The unit is "second". It should be below 1 second. Normally around 0.5s
RR_Tavg = sum(temp)/counter;

% Fetal Heart Rate per minute, physiological FHR should be around 120 to
% 150 beats per minute
FHR_per_min(jj,1) = 60/RR_Tavg;
end

qrs1=[0;qrs1];
FHR_per_min=[FHR_per_min(1);FHR_per_min];
% set(gca,'YColor','white')
% set(gca,'XTick',[0:10:50])
for i=2:length(qrs1)
    line([qrs1(i-1),qrs1(i-1)],[FHR_per_min(i-1),FHR_per_min(i)],'color','b','LineWidth',1.5);
    line([qrs1(i-1),qrs1(i)],[FHR_per_min(i),FHR_per_min(i)],'color','b','LineWidth',1.5);
end
% ylim([120,150]);

% saveFig=1;
% progname=mfilename;
% if(saveFig), figFmt='png';
%         figPath=fullfile('../Figure/',progname);
%         if(~exist(figPath,'dir')), mkdir(figPath); end
%         figName=fullfile(figPath,[cName,'bold_rhythm normal']);
%         print(gcf, ['-d',figFmt],figName);
%         %    saveas(gcf, figName,'fig');
% end

%------------SNR,PNR,CCR等评价指标----------------
ii=2;

start=140;
last=160;
FECG=record(1,(ii-1)*start*fs+1:(ii-1)*start*fs+last*fs)';
ECG=record(2:5,(ii-1)*start*fs+1:(ii-1)*start*fs+last*fs)';
% FECG=record(1,:)';
% ECG=record(2:5,:)';

qrsFcfica1=fix(qrsF);
% qrsFcfica1=fix(qrs);
qrsM1=fix(mother);
qrsM1=qrsM1(qrsM1<=(ii-1)*start*fs+last*fs);
qrsM1=qrsM1(qrsM1>=(ii-1)*start*fs+1);
qrsM1=qrsM1-(ii-1)*start*fs;

qrsFcfica1=qrsFcfica1(qrsFcfica1<=(ii-1)*start*fs+last*fs);
qrsFcfica1=qrsFcfica1(qrsFcfica1>=(ii-1)*start*fs+1);
qrsFcfica1=qrsFcfica1-(ii-1)*start*fs+1;

ref1=FecgDetrFilt(FECG,fs,cName,0,saveFig);
ref1=FecgImpArtCanc(ref1,fs,cName,0,saveFig);
ref=FecgNotchFilt(ref1,fs,cName,0,saveFig);

X=FecgImpArtCanc(ECG,fs,cName,0,saveFig);
Xd=FecgDetrFilt(X,fs,cName,0,saveFig);
Xf=FecgNotchFilt(Xd,fs,cName,0,saveFig);
% Se=FecgICAm(Xf,fs,cName,1,dbFlag,saveFig);



% ------STA估计------------
% spike=zeros(2,start*fs);
% spike(1,qrsM1)=1;
% spike(2,qrsFcfica1)=1;
% cha=1;
% [residual,huifu,huifu_ECG]=bopi(z1,spike,700,1,cha);
% s1=z1(cha,:);
% s2=huifu_ECG(1,:);
% s3=huifu_ECG(2,:);

% ------SVD估计------------
[Xr,Xc]=FecgQRSmCanc(Xf,qrsM1,fs,cName,0,0,saveFig,[]);
% [z1,~,Cond]=yanchibaihua3(Xr',1 ,0);
% Ser1=FecgICAf(z1',fs,cName,graph,dbFlag,saveFig);
% yc=dis_spike(Ser1',0.6);
% [XrF,XcF]=FecgQRSmCanc(yc',qrsFcfica1,fs,cName,0,0,saveFig,[]);
[XrF,XcF]=FecgQRSmCanc(Xr,qrsFcfica1,fs,cName,0,0,saveFig,[]);

for j=1:size(XcF,2)
    if ii==1
        d1=ref(1:150000,1);
        e1=XcF(1:150000,j);
        f1=Xf(1:150000,j);
%         resi=XrF(1:150000,j);
%         m=Xc(1:150000,j);
    else
        d1=ref(10001:160000,1);%refer
        e1=XcF(10001:160000,j);%FECG
        f1=Xf(10001:160000,j);%ORIGION
%         resi=XrF(10001:160000,2);
%         m=Xc(10001:160000,j);
    end
    T=length(d1);

%归一化
d=d1/prctile(abs(d1),99);
e=e1/prctile(abs(e1),99);
f=f1/prctile(abs(e1),99);
% f=f1/16;
%-------------RMSE---------------------
fenzi=min(sum(power(d-e,2)),sum(power(d+e,2)));
RMSE(j)=sqrt(fenzi/T)

fenzi=min(sum(power(d-f,2)),sum(power(d+f,2)));
RMSEyuan(j)=sqrt(fenzi/T)

%--------------correlation coefficient,Pearson相关系数-----------------
% fenzi=sum(d.*e);
% fenmu=sqrt(sum(power(d,2).*sum(e,2)));
temp1=corr(d,e,'type','pearson');
temp2=corr(d,-e,'type','pearson');
CCR(j)=max(temp1,temp2)

temp1=corr(d,f,'type','pearson');
temp2=corr(d,-f,'type','pearson');
CCRyuan(j)=max(temp1,temp2)


%xcorr，是指互相关函数
% temp1=max(xcorr(d,e,'coeff')); 
% temp2=max(xcorr(d,-e,'coeff'));
% Co(j)=max(temp1,temp2)
% 
% temp1=max(xcorr(d,f,'coeff')); 
% temp2=max(xcorr(d,-f,'coeff'));
% Coyuan(j)=max(temp1,temp2)

% %---------------PRD-------------------------
fenzi=min(sum(power(d-e,2)),sum(power(d+e,2)))
fenmu=sum(power(d,2))
PRD(j)=sqrt(fenzi/fenmu)

fenzi=min(sum(power(d-f,2)),sum(power(d+f,2)))
fenmu=sum(power(d,2))
PRDyuan(j)=sqrt(fenzi/fenmu)


%---------------SNR-------------------------
fenzi=sum(power(e,2))
fenmu=min(sum(power(d-e,2)),sum(power(d+e,2)))
SNR(j)=10*log10(fenzi/fenmu)

fenzi=sum(power(f,2))
fenmu=min(sum(power(d-f,2)),sum(power(d+f,2)))
SNRyuan(j)=10*log10(fenzi/fenmu)


end
[~,selectCha]=max(CCR);
% selectCha=1;
Metrics{1,1}=RMSE(selectCha);
Metrics{2,1}=RMSEyuan(selectCha);
Metrics{1,2}=CCR(selectCha);
Metrics{2,2}=CCRyuan(selectCha);
Metrics{1,3}=PRD(selectCha);
Metrics{2,3}=PRDyuan(selectCha);
Metrics{1,4}=SNR(selectCha);
Metrics{2,4}=SNRyuan(selectCha);

% d% ideal ECG,1xN
% e%估计的ECG,1xN
% T=length(d);
% 










