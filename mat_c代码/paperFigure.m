%-------图4 最终结果图-------------
% Power line interference removal by notch filtering 
% Xf=FecgNotchFilt(ECG,fs,cName,0,saveFig);
% [B,A]=butter(2,[5/(fs/2),100/(fs/2)]);
% Xfr=filtfilt(B,A,Xf);%经过filter滤波之后得到的数据y则是经过带通滤波后的信号数据
% Xfr=Xfr(55001:105000,:);
Xfr=Xf(55001:105000,:);

% load('all_r01.mat', 'S')
% % spike=S(:,1:160000);
% spike=S(:,55001:105000);
% [z1,~,Cond]=yanchibaihua3(Xfr',1 ,0);
% residue=bopi(Xfr',spike,700,1,4);

ref1=FecgDetrFilt(FECG,fs,cName,0,saveFig);
ref1=FecgImpArtCanc(ref1,fs,cName,0,saveFig);
ref=FecgNotchFilt(ref1,fs,cName,0,saveFig);

reference=ref(55001:105000);
figure
hold on
set(gcf,'Position',get(0,'ScreenSize'));
set(gcf,'color','w');
set(gca,'Visible','off');
plot(reference,'color',[0.9290 0.6940 0.1250],'linewidth',1);


% load('11-15lastResult2_1.mat');
load('11-15lastResult2_2.mat');

s1=x;
s2=huifu_ECG(1,:);
s3=huifu_ECG(2,:);
resi=x-huifu;
s4=resi;
tm=1/fs:1/fs:50;

% load('11-15lastResult3_1.mat');
load('11-15lastResult3_2.mat');
s5=x;
resi=x-huifu;
s6=resi;


inter=0.03;
ssca=[abs(min(s1))+max(s1);abs(min(s5))+max(s5);abs(min(s2))+max(s2);abs(min(s3))+max(s3);abs(min(s4))+max(s4);abs(min(s6))+max(s6)];
ss=abs(min(s1))+max(s1)+abs(min(s2))+max(s2)+abs(min(s4))+max(s4)+...
    abs(min(s3))+max(s3)+abs(min(s5))+max(s5)+abs(min(s6))+max(s6)+10;
ssca=ssca/ss*(1-inter*6);

nowss=0.01+ssca(1);

%1.recorded AECG signal(channel 1) 

figure
hold on
set(gcf,'Position',get(0,'ScreenSize'));
subplot('Position',[0.01,1-nowss,0.98,ssca(1)]);
nowss=nowss+ssca(2);
hold on
set(gcf,'color','w');
set(gca,'Visible','off');
hold on
plot(tm,-s1,'color',[0 0.4470 0.7410],'linewidth',1.5)
% box off
%-6~6
ylim([min(-s1),max(-s1)])


%2.recorded AECG signal(channel 2) 
% axes(ha(2)); 
subplot('Position',[0.01,1-nowss,0.98,ssca(2)]);
nowss=nowss+ssca(3)+inter-0.01;
hold on
set(gcf,'color','w');
set(gca,'Visible','off');
hold on
plot(tm,s5,'color',[0 0.4470 0.7410],'linewidth',1.5)
% box off
% -11~8.5
ylim([min(s5),max(s5)])


%3.estimianted MECG signal
% axes(ha(3)); 
subplot('Position',[0.01,1-nowss,0.98,ssca(3)]);
nowss=nowss+ssca(4)+inter;
hold on
set(gca,'Visible','off');
plot(tm,-s2,'color',[0.8500 0.3250 0.0980],'linewidth',1.5)
box off
ylim([min(-s2),max(-s2)])
%-1.2~1.5

%4.estimated FECG signal
% axes(ha(4)); 
subplot('Position',[0.01,1-nowss,0.98,ssca(4)+inter/2]);
nowss=nowss+ssca(5)+inter-0.01;
hold on
set(gca,'Visible','off');
hold on
plot(tm,s3,'color',[0.4940 0.1840 0.5560],'linewidth',1.5)
box off
%-0.7~1.4
ylim([min(s3),max(s3)])

%5.residual signal(channel 1) 
% axes(ha(5)); 
subplot('Position',[0.01,1-nowss,0.98,ssca(5)]);
nowss=nowss+ssca(6)+inter-0.01;
hold on
plot(tm,-s4,'color',[0.4660 0.5740 0.1880],'linewidth',1.5)
box off
%-0.6~0.6
ylim([min(-s4),max(-s4)])
set(gca,'Visible','off');


%6.recorded AECG signal(channel 2) 
% -3.6~3.6
% axes(ha(6)); 
subplot('Position',[0.01,1-nowss,0.98,ssca(6)]);
hold on
plot(tm,s6,'color',[0.4660 0.5740 0.1880],'linewidth',1.5)
box off
ylim([min(s6-10),max(s6)])
set(gca,'YColor','white')
% xlabel('Time(s)');
set(gca,'XTick',[0:10:50])

saveFig=1;
% progname=mfilename;
% if(saveFig), figFmt='png';
%         figPath=fullfile('../Figure/',progname);
%         if(~exist(figPath,'dir')), mkdir(figPath); end
%         figName=fullfile(figPath,[cName,'bold_lastResultCon2']);
%         print(gcf, ['-d',figFmt],figName);
%         %    saveas(gcf, figName,'fig');
% end

    
%-------图3 SVD剥离母体心电信号-------------
%debug于FecgQRSmCanc文件142行
%Xx(:,is)-Xc(:,is)
load('figure2Plot.mat');
x1=-Ar(1:400,2);
x1=x1/5;
x2=Xc(55001:70000,2);
x3=Xx(55001:70000,2)-Xc(55001:70000,2);
figure
set(gcf,'color','w');
set(gca,'Visible','off');
tm=1/fs:1/fs:15;
subplot(2,9,[1 2 10 11])
hold on
plot(x1,'linewidth',1.5)
box on
% ylim([min(x1-1),max(x1+1)])
% set(gca,'Visible','off');
ylim([-1,1])
% xlabel('Sample');

subplot(2,9,3:9)
hold on
plot(tm,-x2,'color',[0 0.4470 0.7410],'linewidth',1.5)
box off
ylim([min(x2-1),max(x2)])
% set(gca,'Visible','off');
set(gca,'YColor','white')
% xlabel('Time(s)');

subplot(2,9,12:18)
hold on
plot(tm,x3,'color',[0 0.4470 0.7410],'linewidth',1.5)
box off
ylim([min(x3-1),max(x3)])
set(gca,'YColor','white')
% xlabel('Time(s)');

saveFig=1;
progname=mfilename;
if(saveFig), figFmt='png';
        figPath=fullfile('../Figure/',progname);
        if(~exist(figPath,'dir')), mkdir(figPath); end
        figName=fullfile(figPath,[cName,'bold_SVDCanc']);
        print(gcf, ['-d',figFmt],figName);
        %    saveas(gcf, figName,'fig');
end
    


%-------图2 识别母体心电信号-------------

yuan=ECG(35001:50000,3);
yuan2=ECG(35001:50000,4);
yuchuli=-Xf(35001:50000,3);
yuchuli2=Xf(35001:50000,4);
fastica=-Se(35001:50000,1);
tm=1/fs:1/fs:15;

s1=ECG(35001:50000,3);
s2=ECG(35001:50000,4);
% s3=-Se(35001:50000,1);
% s4=ECG(35001:50000,2);
s3=-Xf(35001:50000,3);
s4=Xf(35001:50000,4);
s5=-Se(35001:50000,1);


inter=0.03;
ssca=[abs(min(s1))+max(s1);abs(min(s2))+max(s2);abs(min(s3))+max(s3);abs(min(s3))+max(s4);abs(min(s4))+max(s3);abs(min(s3))+max(s3)];
ss=abs(min(s1))+max(s1)+abs(min(s2)+max(s2))+abs(min(s3)+max(s3))+(abs(min(s4))+max(s4))*3;
ssca=ssca/ss*(1-inter*6);
nowss=0.01+ssca(1);

figure
hold on
set(gcf,'Position',get(0,'ScreenSize'));

%原信号1
% subplot(8,1,[1 2])
subplot('Position',[0.01,1-nowss,0.98,ssca(1)]);
hold on
nowss=nowss+ssca(2);
set(gcf,'color','w');
set(gca,'Visible','off');
plot(tm,s1,'color',[0 0.4470 0.7410],'linewidth',1.5)
box off
ylim([min(s1),max(s1)])

%原信号2
% subplot(8,1,[1 2])
subplot('Position',[0.01,1-nowss,0.98,ssca(2)]);
hold on
nowss=nowss+ssca(3)+inter/2;
set(gcf,'color','w');
set(gca,'Visible','off');
plot(tm,s2,'color',[0 0.4470 0.7410],'linewidth',1.5)
box off
ylim([min(s2),max(s2)])

%预处理后的信号
% subplot(8,1,[3 4])
subplot('Position',[0.01,1-nowss,0.98,ssca(3)]);
nowss=nowss+ssca(4)+inter/2;
set(gca,'Visible','off');
hold on
plot(tm,-s3,'color',[0.8500 0.3250 0.0980],'linewidth',1.5)
box off
ylim([min(-s3),max(-s3)])

%预处理后的信号2
% subplot(8,1,[3 4])
subplot('Position',[0.01,1-nowss,0.98,ssca(4)]);
nowss=nowss+ssca(5)+inter/2;
set(gca,'Visible','off');
hold on
plot(tm,s4,'color',[0.8500 0.3250 0.0980],'linewidth',1.5)
box off
ylim([min(s4),max(s4)])


%FastICA后的信号
% subplot(8,1,[5 6])
subplot('Position',[0.01,1-nowss,0.98,ssca(5)]);
nowss=nowss+ssca(6)+inter/2;
set(gca,'Visible','off');
hold on
plot(tm,fastica,'color',[0.4940 0.1840 0.5560],'linewidth',1.5)
box off
ylim([min(fastica),max(fastica)])
line([1/fs,15],[4,4],'linestyle','--','color','r');


%detected peaks
% subplot(8,1,[7 8])
subplot('Position',[0.01,1-nowss,0.98,ssca(6)]);
hold on
% plot(tm,fastica,'color',[0.4660 0.6740 0.1880])
box off
ylim([-0.5,1])
[~,L2] = findpeaks(fastica,'minpeakheight',4,'MINPEAKDISTANCE',200);
for i=1:length(L2)
    line([L2(i)/fs,L2(i)/fs],[0,1],'linestyle','-','LineWidth',2,'color',[0.6350 0.0780 0.1840]);
end
set(gca,'YColor','white')
% xlabel('Time(s)');

saveFig=1;
progname=mfilename;
if(saveFig), figFmt='png';
        figPath=fullfile('../Figure/',progname);
        if(~exist(figPath,'dir')), mkdir(figPath); end
        figName=fullfile(figPath,[cName,'bold_MECG_CFICA']);
        print(gcf, ['-d',figFmt],figName);
        %    saveas(gcf, figName,'fig');
end
    

%-------图1 cfica和ica的对比-------------
load('paperFigure1.mat');

fs=2000;
figure
set(gcf,'Position',get(0,'ScreenSize'));
set(gcf,'color','w');
set(gca,'Visible','off');
s1=yc(2,80001:110000);
s2=yt(1,80001:110000);
tm=1/fs:1/fs:15;

% FECG参考信号
saveFig=0;
ref1=FecgImpArtCanc(FECG,fs,cName,0,saveFig);
ref1=FecgDetrFilt(ref1,fs,cName,0,saveFig);
ref1=FecgNotchFilt(ref1,fs,cName,0,saveFig);
ref=ref1(40001:55000,1)';
tm1=1/1000:1/1000:15;

subplot(14,1,[1 3])
set(gca,'Visible','off');
hold on
plot(tm1,ref,'color',[0.6350 0.0780 0.1840],'linewidth',1.5)
box off
ylim([min(ref),max(ref)])

subplot(14,1,[4 6])
set(gca,'Visible','off');
hold on
plot(tm,s1,'linewidth',1.5)
box off
ylim([min(s1),max(s1)])
line([1/fs,15],[2.3,2.3],'linestyle','--','color','r');

subplot(14,1,[7 8])
hold on
[d1,L1] = findpeaks(s1,'minpeakheight',2.3,'MINPEAKDISTANCE',200);
for i=1:length(d1)
    line([L1(i)/fs,L1(i)/fs],[-1,1],'linestyle','-','LineWidth',2,'color',[0.4940 0.1840 0.5560]);
end
set(gca,'Visible','off');
box off

subplot(14,1,[9 11])
set(gca,'Visible','off');
hold on
plot(tm,s2,'color',[0.4660 0.5740 0.1880],'linewidth',1.5)
box off
ylim([[min(s2),max(s2)]])
line([1/fs,15],[3.0,3.0],'linestyle','--','color','r');


subplot(14,1,[12 13])
hold on
[d2,L2] = findpeaks(s2,'minpeakheight',3.1,'MINPEAKDISTANCE',200);
for i=1:length(d2)
    line([L2(i)/fs,L2(i)/fs],[-1,1],'linestyle','-','LineWidth',2,'color',[0.8500 0.3250 0.0980]);
end
set(gca,'Visible','off');
box off

subplot(14,1,14)
hold on
box off
plot(tm,1)
set(gca,'YColor','white')
% xlabel('Time/s');

saveFig=1;
progname=mfilename;
if(saveFig), figFmt='png';
        figPath=fullfile('../Figure/',progname);
        if(~exist(figPath,'dir')), mkdir(figPath); end
        figName=fullfile(figPath,[cName,'bold_FECG_CFICA']);
        print(gcf, ['-d',figFmt],figName);
        %    saveas(gcf, figName,'fig');
end
