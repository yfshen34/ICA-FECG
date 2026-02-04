clc
clearvars -except Qm Qf wm wf
% clear all
% close all
inter=2;
%  for user = 5
%     for user =[3,4,5,7,14,20,27,34,40,47,51,58,63,68,73]
user=5;
%  try

if user<10
    fname=strcat('a0', num2str(user));
    name=strcat('a0', num2str(user),'.fqrs.txt');
else
    fname=strcat('a', num2str(user));
    name=strcat('a', num2str(user),'.fqrs.txt');
end

[signal,fs]= rdsamp(fname,[], []);

ann=rdann(fname,'fqrs',[],[]);
ECG0=signal(10001:end,:);
qrsAf=ann;
% load('b1515.mat');
% signal=data15(20001:56000,:);
% ECG=signal;fs=200;

t=[];
Xm0=[];Xf0=[];

XR=[];
Xrf=[];
qrsM_all=[];
qrsF_all=[];
fhr_seg = [];
step = 50;
step_on = 1001;
step_off = 1050;
window =3000;

imax = floor((50000 - window)/step);
m=1;

for i=0:imax
    ECG1=ECG0(step*i+1:step*i+window,:);


    ta=cputime;
    tm=0.001:0.001:60;
    tm=tm';
   
    
%     cName=fname;
%    
%     qrsAf=ann;
   % t1=(o-19)/fs:1/fs:o/fs;
  %  t=[t t1];
  cName=0;

    dbFlag=0;                   % debug flag
    graph=0;                    % enable/disable graphical representation
    saveFig=0;                  % =1 => save figures of the processing phases
    saveFigRRf=0;               % =1 => save estimated fetal RR figures
    
    X=FecgImpArtCanc(ECG1,fs,cName,0,saveFig);
    % close all
    
%     ---- detrending  ----低通滤波
    Xd=FecgDetrFilt(X,fs,cName,graph,dbFlag);
%     Xd=FecgDetrFilt(X,fs,cName,0,saveFig);
%     close all
    
%     ---- Power line interference removal by notch filtering ----
    Xf=FecgNotchFilt(Xd,fs,cName,graph,dbFlag);
%     Xf=FecgNotchFilt(Xd,fs,cName,0,saveFig);
    
%     Se=Qm*Xf';
    Se=Qm*X';
    Se=wm*Se;
    Se=Se';

%     Xr1=Xf';
     Xr1=X';
%     [Nsignal, Nsample] = size(Xf');
    [Nsignal, Nsample] = size(X');
    R=1;
    xt = zeros((2*R+1)*Nsignal,Nsample);
    for k= 1:Nsignal
        %xt(R*(i-1)+1:R*i,:)=toeplitz([x(i,1);zeros(R-1,1)],x(i,:)); toeplitz(x)用向量x生成一个对称的托普利兹矩阵
        xt((2*R+1)*(k-1)+1:(2*R+1)*k,:)=toeplitz([Xr1(k,(R+1):-1:1),zeros(1,R)],[Xr1(k,(R+1):end),zeros(1,R)]);%xt为扩展矩阵
    end
    xt=xt-mean(xt,2)*ones(1,Nsample);%将xt中心化 ,求每行的均值

  
 
    % ---- Channel selection and Mother QRS detection
    [qrsM,indexm]=FecgQRSmDet(Se,fs,cName,0,dbFlag,saveFig,qrsAf);


    if i == 0
    qrsM_step_loc = find(qrsM <= step_off);
    elseif i == imax
    qrsM_step_loc = find(qrsM >= step_on);
    else
    qrsM_step_loc = find(qrsM <=step_off & qrsM >= step_on);
    end
    qrsM_all = [qrsM_all; qrsM(qrsM_step_loc) + i * step];



    % ---- Mother QRS cancelling
    [Xr,Xm]=FecgQRSmCanc(X,qrsM,fs,cName,1,1,saveFig,qrsAf);
    
    if(i==0)
        XR=Xr;
    end
    % XR=[XR(51:1000,:);Xr(951:1000,:)];
    % Xm=X-XR;
    Xm1=Xm(step_on:step_off,:);
    Xm0=[Xm0;Xm1];

    % Xm1=Xm(:,indexm);
%      for vi=1:size(Xm,2)
%         a=diff(Xm(:,vi));
%         % much=(length(a(a<1000))>=sum(a(a<1000))/fs*3);%发放个数的限制
%         a((a>500))=[];
%         b=Xm(qrsM,vi);
%         coef_s1(vi)=std(a)/mean(a);
%         coef_w1(vi)=std(b)/abs(mean(b));
%     end
%     [cs1,idsm]=sort(coef_s1);%%s interval
%       [cw1,idw]=sort(coef_w1);%% w amp
%     sw=coef_s1+0.3*coef_w1;
%     [csw,idsw]=sort(sw);
%     Xm1=Xm(:,idsw(1));

%     figure;
%     % plot(Xm(:,2));
%      stackedplot(Xm);
% % 
% %     % hold on;
% %     % plot(960,Xm(960,1),'.','Color','r','MarkerSize',25);
%       pause(1);
% 
%      Xr1=derivative(Xr,fs);
%      % Xr1=Xr1';
%      Xr1=Xr';
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

    
    [qrsFcfica,indexf]=FecgQRSfDet(Ser1,fs,cName,qrsM,graph,dbFlag,0,saveFigRRf,qrsAf);
    
    qrsf_win=qrsFcfica*fs;


if i == 0
    qrsf_step_loc = find(qrsf_win <= step_off);
elseif i == imax
    qrsf_step_loc = find(qrsf_win >= step_on);
else
    
    qrsf_step_loc = find(qrsf_win <=step_off & qrsf_win >= step_on);
end

qrsF_all = [qrsF_all; qrsf_win(qrsf_step_loc) + i * step];

%     [M,N]=size(qrsf_win);
%     RRf_seg(m,1)=(qrsf_win(M,1)-qrsf_win(1,1))/(M-1)/fs;
%     m=m+1;
% 
%     RR_f = diff(qrsf_win)/fs;
%     fhr=60/mean(RR_f);
%     fhr_seg = [fhr_seg;fhr];

    [~,Xf]=FecgQRSfCanc(Xr,qrsf_win,fs,cName,1,1,saveFig,qrsAf);
    Xf1=Xf(step_on:step_off,:);
    Xf0=[Xf0;Xf1];
    Xrf=[Xrf;Ser1(501:550,:)];
  
%   stackedplot(Xf);
%   pause(1);
%     close all;

    tb=cputime-ta;
     disp('online运行时间: ');
     disp(tb);
     close all;
%     % Xf=Xf(:,indexf);
%      for i=1:size(Xf,2)
%         a=diff(Xf(:,i));
%         % much=(length(a(a<1000))>=sum(a(a<1000))/fs*3);%发放个数的限制
%         a((a>500))=[];
%         b=Xf(floor(qrsF),i);
%         coef_s2(i)=std(a)/mean(a);
%         coef_w2(i)=std(b)/abs(mean(b));
%     end
%     [cs1,idsf]=sort(coef_s2);%%s interval
%       [cw1,idw]=sort(coef_w2);%% w amp
%     sw=coef_s2+0.3*coef_w2;
%     [csw,idsw]=sort(sw);
%     Xf1=Xf(:,idsw(1));
%     figure;plot(Xf1);
%     close all;
 end
% err=ans-qrsAf(44:65,1);
% mean(abs(err));
% RR_qrsAf=diff(qrsAf)/1000;
% RR_qrsF_all=diff(qrsF_all)/fs;

qrsAf=qrsAf/1000;
qrsF_all=(qrsF_all)/fs;
% 
% 
% 
% 
%  indexf_all(i)=indexf;
%  qrsFcfica_all(i)=qrsFcfica;
% 
% 
%  filename3=['sub',num2str(user),'_qrsFcfica.txt'];
% % 
% %  target_dir = 'E:\mat_c代码 _5.11\mat_c代码\fDet_test_data';
% %  if ~exist(target_dir, 'dir')
% %     mkdir(target_dir);  % 创建目录
% %     if ~exist(target_dir, 'dir')  % 再次检查是否创建成功
% %         error('无法创建目录: %s', target_dir);  % 如果仍失败，报错退出
% %     end
% % end
% 
% filename1=['fDet_test_data\sub',num2str(user),'_Ser1.txt'];
% filename2=['fDet_test_data\sub',num2str(user),'_qrsM.txt'];
% 
% Filename1=['fDet_test_data\sub',num2str(user),'.mat'];
% % % save(Filename1, 'Serx','qrsM_all'); 
% 
% % % save(filename1, 'Serx','-ascii', '-tabs');  % 使用 -tabs 指定制表符分隔
% save(filename2, 'qrsM_all', '-ascii',"-tabs");  
% 
% 
% 
% % save(filepath1, 'Ser1', '-ascii', '-tabs');  % 使用 -tabs 指定制表符分隔
% % save(filepath2, 'qrsM_all', '-ascii',"-tabs");  
% % save(filepath3, 'qrsFcfica_all', '-ascii',","); 
% 
% catch ME
%     % 如果发生错误，打印错误信息并继续下一个循环
%     warning('Error processing user %d (%s): %s', user, fname, ME.message);
%     continue;  % 跳过当前用户，继续下一个循环
% end
% end



