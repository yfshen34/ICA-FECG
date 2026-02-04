
clc
% clear
clearvars -except Qm Qf wm wf
close all
t1=cputime;
load("不同水平仿真数据.mat");
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
ann=floor(ann/5);
ECG=signal;
tm=0.001:0.001:60;
tm=tm';
cName=fname;
% qrsAf=ann;
qrsAf=[];

dbFlag=0;                   % debug flag
graph=0;                    % enable/disable graphical representation
saveFig=0;                  % =1 => save figures of the processing phases
saveFigRRf=0;               % =1 => save estimated fetal RR figures

fs = 1000;             % sampling frequency
ECG=load("G:/data/05-11-1515-filter.txt");
% load("testdata1.mat");
% ECG=emgdata;
ECG=ECG';

% [b,a0] = butter(5,[20/(fs/2),250/(fs/2)]); 
% % 
% ECG=filtfilt(b,a0,ECG);
% 
% a2 = [1,-3.05311331771918e-16,2.62355180660524,-6.08122226439088e-16,2.31468258108911,-2.86629396097301e-16,0.685535977284663];
% b2 = [0.827971295622376,-3.04191719090949e-16,2.48391388686713,-6.08383438181898e-16,2.48391388686713,-3.04191719090949e-16,0.827971295622376];
% ECG=filtfilt(b2,a2,ECG);
% ECG=load("G:/data/329cf1.txt");
% ECG=test5;
% ECG=double(xinhao5');
ECG=ECG(1:12000,:);
% ECG=resample(ECG,1,5);
fs=200;
 
% ECG=load("data.txt");
% ECG=test7;

% ECG=ECG(1:12000,:);

% X=FecgFecgImpArtCanc(ECG,fs,cName,graph,dbFlag);
% X=FecgDetrFilt(ECG,fs,cName,0,saveFig);
% a1 = [1,-0.0399789305697524,-1.95035734571979,0.0380715113700012,0.952303292302367];
% b1 = [0.975860239678118,0,-1.95172047935624,0,0.975860239678118];
% a2 = [1,-3.05311331771918e-16,2.62355180660524,-6.08122226439088e-16,2.31468258108911,-2.86629396097301e-16,0.685535977284663]
% b2 = [0.827971295622376,-3.04191719090949e-16,2.48391388686713,-6.08383438181898e-16,2.48391388686713,-3.04191719090949e-16,0.827971295622376]
% X1=filtfilt(b1,a1,ECG);
% X2=filtfilt(b2,a2,X1);





% ---- Artifact canceling ----
% 去尖峰噪声
X=FecgImpArtCanc(ECG,fs,cName,0,saveFig);
% X=ECG;
% 
% a1 = [1,-0.0399789305697524,-1.95035734571979,0.0380715113700012,0.952303292302367];
% b1 = [0.975860239678118,0,-1.95172047935624,0,0.975860239678118];
% [b,a0] = butter(1,[30/(fs/2),99/(fs/2)]); 
% 
% X=filtfilt(b1,a1,X);


% close all
% 
% while(max(X,[],"all")<=50 || max(X,[],"all")>=200)
%     if max(X)<=50
%         X = X*2;
%     elseif max(X)>=200
%         X = X/2;
%     end
% end
% X=X*100;

Xr = X;
qrsM=0;
Qm=zeros(48,48);
Qf=zeros(48,48);
wm=zeros(4,4);
wf=zeros(4,4);
for k=1:3
        
    X0=Xr;
    po=0;
    % if k==1
    %     po=0.5;
    % else
    %     po=0;
    % end
    [Xf,Q,Cond1]=yanchibaihua3(X0',1 ,po);
    
    
    
   % ---- Independent Component Analysis ----
   % Xm=FecgICAm(Xf,fs,cName,graph,dbFlag,saveFig);
    Se=FecgICAm(Xf',fs,cName,0,dbFlag,saveFig);

    yc=dis_spike(Se',0.6);

    peakinterval=0.3*fs;
    debug=0;
    % ann1=ann/1000*fs;
    ann1=[];
%     yc=[Se yc']';
    tmp=yc';
    F=[];
    coef_s=[];
    coef_w=[];
    ind=[];
    for i=1:size(yc,1)
        [F1,yc(i,:),th,coef_s1,coef_w1,ind1]=getspike(yc(i,:),fs,peakinterval,debug,ann1);
        coef_s=[coef_s,coef_s1];
        coef_w=[coef_w,coef_w1];
        ind=[ind,ind1];
        F=[F,F1];
    end
     [cs1,ids]=sort(coef_s);%%s interval
    [cw1,idw]=sort(coef_w);%% w amp
    sw=coef_s+0.3*coef_w;
    [csw,idsw]=sort(sw);
    for i=1:min(4,size(yc,1))

        a{1,i}=F{1,idsw(i)}; 

    end
    F=a;
    
    % for i=1:size(F,2)
    %     a=diff(F{i});
    %     % much=(length(a(a<1000))>=sum(a(a<1000))/fs*3);%发放个数的限制
    %     a((a>500))=[];
    %     b=yc(i,F{i});
    %     coef_s(i)=std(a)/mean(a);
    %     coef_w(i)=std(b)/abs(mean(b));
    % end
    % F_new = cell(1,size(F,2));
    % F_newdiff = [];
    % F_new = [];
    % Fm = [];
    % for i=1:size(F,2)
    %     % Fdiff1 = diff(F{i});
    %     % F_new1  = Fdiff1(1:50);
    %     % F_newdiff   = [F_newdiff;F_new1];
    %     % F_new2  = F{i}(1:50);
    %     % F_new   = [F_new;F_new2];
    %     % F_new3  = mean(Fdiff1);
    %     % Fm(i,1) = F_new3;
    %     Fm(i,2) = size(F{i},2);
    %     % Fdiff2 = diff(Fdiff1);
    %     % 
    %     % F_index = sum(abs(Fdiff2));
    %     % Fm(i,1) = sum(abs(Fdiff1));
    % 
    %     % Ftest(i) = F_index;
    % 
    %     % plot(Fdiff1, 'LineWidth', 2,'DisplayName',['Fdiff1_=',num2str(i),',index=',num2str(F_index)]);
    %     % hold on; grid on;
    %     % pause(0.3);
    %     % legend('Location', 'best', 'AutoUpdate', 'on');
    % 
    % end
    % % [F_A,F_I] = sort(-Ftest);
    % 
    % 
    % [num,eq] = hist(Fm(:,2),5);
    % [num_num,num_ind] = sort(-num);
    % width = eq(2) - eq(1);
    % Fleft  = eq(num_ind(1)) - width/2;
    % Fright = Fleft + width;
    % j=1;
    % for i = 1:size(F,2)
    %     if (size(F{i},2)>=Fleft && size(F{i},2)<=Fright)
    %         Fdiff1 = diff(F{i});
    %         Fdiff2 = diff(Fdiff1);
    %         F_index = sum(abs(Fdiff2));
    %         Ftest(j) = F_index;
    %         F_I(j) = i;
    %         j = j+1;
    %     end
    % end
    % [F_A,F_I_index] = sort(Ftest);
    % final = min(2,size(F,2));
    % F_final = cell(1,final);
    % for i=1:final
    %     F_final{i} = F{F_I(F_I_index(i))};
    % end
    % F = F_final;
    % a{1} = F{3};
    % F = a;
    win=0.01*fs;
    if ~isempty(F)
        refer=zeros(size(F,2),size(Se,1));
        for i=1:size(F,2)
            if ~isempty(F{i})
                for j=1:size(F{i},2)
                    refer(i,max(F{i}(j)-win,1):min(F{i}(j)+win,size(Se,1)))=1;
                end
            end
            
        end
    end
    i=1;
    yt=[];
    w=[];
    while i<=size(refer,1)
         [ytt, wtt] = cfICA(Xf, refer(i,:),qrsM,qrsAf,fs,user);
         yt(:,i)=ytt';
         w(:,i)=wtt;
         i=i+1;
    end
    
    if(k<=1)
        [qrsM,indexm]=FecgQRSmDet(yt,fs,cName,0,dbFlag,saveFig,qrsAf);
        [Xr,Xm]=FecgQRSmCanc(X0,qrsM,fs,cName,0,0,saveFig,qrsAf);%%Xm为母心信号      
        Qm=Q;
        wm=w';
    else       
        qrsFcfica=FecgQRSfDet(yt,fs,cName,qrsM,0,0,0,saveFigRRf,qrsAf);
        qrsf{k-1}=floor(qrsFcfica*fs);
        Qf=Q;
        wf=w';
        [Xr,Xfee]=FecgQRSfCanc(X0,qrsf{k-1},fs,cName,0,0,saveFig,qrsAf);
        Xfe{k-1}=Xfee;%%Xfe为胎心信号
    end
end
t2=cputime-t1;
disp('offline运行时间: ');
disp(t2);
