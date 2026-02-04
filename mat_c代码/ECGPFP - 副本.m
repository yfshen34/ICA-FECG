function [SpikeC,W,residue,huifu_ECG,huifu,peelcount,numc]=ECGPFP(x,Spike,param,ann)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明

% ---- paramater setting ----
if nargin==0
end
if nargin==1
    Spike=[];
    param=struct();
end
if nargin==2
    param=struct();
end
if (~isfield(param,'peelnum')) % isfield判断输入是否是结构体数组的域（成员）。
    param.peelnum=3;
end
if (~isfield(param,'ficanum'))
    param.ficanum=5;
end
if (~isfield(param,'outnum'))
    param.outnum=5;
end
if (~isfield(param,'fs'))
    param.fs=2048;
end
if (~isfield(param,'wavelength'))
    param.wavelength=fix(param.fs/20);
end
if (~isfield(param,'peakinterval'))
    param.peakinterval=60/200*param.fs;
end

if (~isfield(param,'df'))
    param.df=5;
end
if (~isfield(param,'dc'))
    param.dc=20;
end
if (~isfield(param,'CORR'))
    param.CORR=0.5;
end
if (~isfield(param,'covs'))
    param.covs=0.4;
end
if (~isfield(param,'covp'))
    param.covp=0.35;
end
if (~isfield(param,'fig'))
    param.fig=1;
end
if (~isfield(param,'auto'))
    param.auto=1;
end

if (~isfield(param,'ini')) %
    param.ini=0;
end

if (~isfield(param,'out')) %
    param.out=100;
end

coef=0.6;

% ---- residual  ----
x=x-mean(x,2)*ones(1,size(x,2));%x中每个元素减去均值
x=single(x);%signal四个字节，节省空间。
if isempty(Spike) %判断输入是否为空
    residue=x;
else
    residue=bopi(x,Spike,param.wavelength,0);
end

% ---- parameter ----
count=1;

SpikeC=[];
W=[];
huifu_ECG=[];
huifu=[];

% yc=zeros(param.ficanum*param.outnum,size(x,2));%剥皮前FastICA运行数*输出数，信号采样数
yc=[];
%tongdao=randperm(size(x,1));
%xw=yanchibaihua3(x(tongdao(1:min(size(x,1),64)),:),param.dc,0.0);
% xw=yanchibaihua3(x,param.dc,0.0);

% ---- starting  ----
while count<=param.peelnum
    % ---- delayed whiten  ----
    ran=randi(param.df)
    [z,~,Cond]=yanchibaihua3(residue,param.df ,0);
    if Cond>1e15
        disp(['warn:the Cond is too large!',num2str(Cond)]);
        return
    end
    
    % ----deflation FastICA ----
    for i=1:param.ficanum
        [Se,Ae,cerr]=coshFpDeIca(z(1:floor((0.75+0.25*(rand>0.75))*size(z,1)),:),param.outnum);
        yc=[yc;Se];
%         yc(((i-1)*param.outnum+1):i*param.outnum,:)=Se;
        
        %         if param.debug
        %         multiPlot('FastICA',Se);
        %         end
        
    end
    % delete same yc
    yc=dis_spike(yc,coef);
    
    
    % ----threshold selecting ,clustering----
    F=[];
    for i=1:size(yc,1)
        [F1,F2,yc(i,:),th]=getspike(yc(i,:),param.fs,param.peakinterval,param.debug);
        F=[F,F1];
    end
    
    % generate reference signal,20ms
    win=0.01*param.fs;
    if ~isempty(F)
        refer=zeros(size(F,2),size(x,2));
        for i=1:size(F,2)
            if ~isempty(F{i})
                for j=1:size(F{i},2)
                    refer(i,F{i}(j)-win:F{i}(j)+win)=1;
                end
            end
            
        end
    end
    
    
    % ----constraint fastICA ----
    xw=yanchibaihua3(x,param.dc,0.0);
    i=1;
    while i<=size(refer,1)
%         yc(i,:)=fcica3(xw,refer(i,:));
        [yc(i,:), ~] = cfICA(xw, refer(i,:));
        ycmax = max(yc(i,:));
        figure;set(gcf,'Position',get(0,'ScreenSize'));
        plot(yc(i,:));title([num2str(i),'/',num2str(size(yc,1)),'Please input the threshold(or press "enter" to auto setting or input "0" to delete it):']);
        
        %——————输入————
        while(1)
            de=input('Please input the threshold(or press "enter" to auto setting or input "0" to delete it):','s');
            if strcmp(de ,'i')||strcmp(de ,'j')
                continue;
            end
            if (isempty(de))
                break;
            else
                de = str2num(de);
                if (~isempty(de))&& (de<ycmax)
                    break;
                end
            end
        end
        
        close
        
        if de==0 %输入0，结束这次循环
            i=i+1;
            continue;
        end
        
        %————为空自动选择阈值，否则按输入阈值————
        if isempty(de)
            TH=spikethresh(yc(i,:));
            S=searchspike_in(yc(i,:),TH,param.peakinterval);
%             yt=fcica3(xw,S);
            [yt, ~] = cfICA(xw, S);
        else
            S=searchspike_in(yc(i,:),de,param.peakinterval);
%             yt=fcica3(xw,S);
            [yt, ~] = cfICA(xw, S);
        end
        ytmax = max(yt);
        figure;set(gcf,'Position',get(0,'ScreenSize'));
        subplot(211);plot(yt);title('Determine the threshold (or "enter" to accept or "0" to discard):');
        subplot(212);plot(S);
        hold on
        plot(ann,1,'r+');

        
        %————第一次输入阈值，接受/再次输入阈值————
        while(1)
            de=input('Determine the threshold (or "enter" to accept or "0" to discard):','s');
            if strcmp(de ,'i')||strcmp(de ,'j')
                continue;
            end
            %————接受————
            if (isempty(de))
                break;
            else
                %————重新输入阈值————
                de = str2num(de);
                if (~isempty(de))&& (de<ytmax)
                    break;
                end
            end
        end
        close
        
        %————重新输入的是数字————
        while (~isempty(de) & de~=0)
            S=searchspike_in(yt,de,param.peakinterval);
%             yt=fcica3(xw,S);
            [yt, ~] = cfICA(xw, S);
            ytmax = max(yt);
            figure;set(gcf,'Position',get(0,'ScreenSize'));
            subplot(211);plot(yt);title('Determine the threshold (or "enter" to accept or "0" to discard):');
            subplot(212);plot(S);
            hold on
        plot(ann,1,'r+');
            
            %————接受/重新设置阈值————
            while(1)
                de=input('Determine the threshold (or "enter" to accept or "0" to discard):','s');
                if strcmp(de ,'i')||strcmp(de ,'j')
                    continue;
                end
                %————接受————
                if (isempty(de))
                    break;
                else
                    %————否则输入一个数字————
                    de = str2num(de);
                    if (~isempty(de))&& (de<ytmax)
                        break;
                    end
                end
            end
            close
            
        end
        
        %接受Spike，进行剥离(母体心电)，或者胎儿心电估计(STA)
        if de==0
            S=[];
        else
            
            %————peel off————
            if isempty(Spike)
                % SVD剥去MECG
                residue=FecgQRSmCanc(x,S,param.fs,[],param.debug,0,0,ann);
                xw=FecgQRSmCanc(xw,S,param.fs,[],param.debug,0,0,ann);
                
                %————加入序列————
                Spike=[Spike;S];
                
            elseif ~cmp_spike(Spike,S,0.5)
                % %判断与之前的spike是否重复不重复，STA 波形估计MECG,FECG
                Spike=[SpikeC;Spike];
                [residue,huifu,huifu_ECG,W]=bopi(x,Spike,param.wavelength,1);
                
                %————加入序列————
                Spike=[Spike;S];
                
                %————分解结束————
                SpikeC=Spike;
                numc=size(SpikeC,1);
                peelcount=count;
                return;
            else
                fprintf('重复Spike，无效\n');
            end
            
            
        end        
        
    end
    %————是否继续分解————
            title('Press "enter" to continue or "0" to finish):')
            re=input('Press "enter" to continue or "0" to finish):');
            if ~isempty(re)||count==param.peelnum
                if ~isempty(Spike)
                    [residue,huifu,huifu_ECG,W]=bopi(x,Spike,param.wavelength,1);
                
                SpikeC=Spike;
                numc=size(SpikeC,1);
                peelcount=count;
                return;
                end
            end
    % ————比较是否出现重复序列————
%     Spike=dis_spike(Spike,0.2);
%     peelcount(count)=size(Spike,1);
%     %若没有spike则开始下一次FastICA
%     if isempty(Spike)
%         count=count+1;
%         continue;
%     end
    
    % ————排序————
%     temp=[];
%     for j=1:size(Spike,1)
%         temp(j)=length(find(Spike(j,:)==1));
%     end
%     if size(Spike,1)>=2
%         if temp(1)>temp(2)
%             temp_S=Spike(1,:);
%             Spike(1,:)=Spike(2,:);
%             Spike(2,:)=temp_S;
%         end
%     end
    
    %     % ----peel off----
    % %         [residue,huifu,huifu_ECG,W]=bopi(x,Spike,param.wavelength,1);
    %     %————现有的Spike数量————
    %     if size(Spike,1)==1
    %         if isempty(SpikeC)
    %             % SVD剥去MECG
    %                 [residue,huifu,huifu_ECG,W]=FecgQRSmCanc(x,Spike,param.fs,[],param.debug,0,0,ann);
    %                 SpikeC=[SpikeC;Spike];
    %         else
    %             % STA 波形估计MECG,FECG
    %         Spike=[SpikeC;Spike];
    %         [residue,huifu,huifu_ECG,W]=bopi(x,Spike,param.wavelength,1);
    %         end
    %     else
    %         % STA 波形估计MECG,FECG
    %         [residue,huifu,huifu_ECG,W]=bopi(x,Spike,param.wavelength,1);
    %     end
%     title('Press "enter" to continue or "0" to finish):')
%     re=input('Press "enter" to continue or "0" to finish):');
%     %不为空，继续下一次分解
%     if ~isempty(re)
%         SpikeC=Spike;numc=size(SpikeC,1);
%         return;
%     end
%     %达到剥皮次数上限，退出
%     if count==param.peelnum
%         return;
%     else
%         %             close
%     end
%      
    count=count+1;
    
end


    function [S]=dis_spike(y,p,fs,graph)
        if isempty(p)
            p=0.4;
        end
        if nargin<3
            fs=1000;
        end
        if nargin<4
            graph=param.debug;
        end
        S=[];
        win=fs/200;
        while ~ isempty(y) %isempty判断输入是否为空
            S=[S;y(1,:)];
            Co=zeros(size(y,1),1);
            for j=1:size(y,1)
                temp1=max(xcorr(y(1,2*win+1:end-2*win),y(j,2*win+1:end-2*win),'coeff')); %xcorr，是指互相关函数
                temp2=max(xcorr(y(1,2*win+1:end-2*win),-y(j,2*win+1:end-2*win),'coeff'));
                Co(j)=max(temp1,temp2);
            end
            y(Co>=p,:)=[];
        end
        
        if graph
            multiPlot('shaixuanguohoudey',S);
        end
        
    end

%寻找当前spike是否为重复变量
    function [flag]=cmp_spike(y,s,p)
        if nargin<3
            p=0.2;
        end
        
        win=param.fs/200;
        for j=1:size(y,1)
            temp1=max(xcorr(s(1,2*win+1:end-2*win),y(j,2*win+1:end-2*win),'coeff')); %xcorr，是指互相关函数
            temp2=max(xcorr(s(1,2*win+1:end-2*win),-y(j,2*win+1:end-2*win),'coeff'));
            Co(j)=max(temp1,temp2);
        end
        if find(Co>=p)
            flag=0;
        else
            flag=1;
        end
        
    end

%————找到x中峰值的位置————
    function S=searchspike_in(x,high,peak_interval,fs)
        if nargin<4
            fs=param.fs;
        end
        
        if nargin<5
            graph=param.debug;
        end
        
        S=zeros(size(x));
        win1=fix(fs/200);
        x_jd=x(2*win1+1:end-2*win1);
        
        % 去趋势
        in=x_jd';
        fin=[in(1:fs);in];
        np = length(fin);
        xd = decimate(fin,round(fs/50),'fir'); % decimation
        lbx = medfilt1(xd,10); % median filtering - get isoline
        lb = interp(lbx,round(fs/50)); % interpolation
        out = fin-lb(1:np);
        out(1:fs)=[];  % result
        X=out';
        
        if graph
            figure
            subplot(2,1,1); hold on;
            plot(x_jd)
            legend('选取constrainted FastICA的spike去趋势前：');
            plot(lb(fs+1:np),'r')
            legend('基线');
            
            subplot(2,1,2),hold on
            plot(X)
            legend('选取constrainted FastICA的spike去趋势后：');
            hold off
        end
        
        if high==0
            S=[];
            return
        elseif high>0
            yp=x_jd;
        else
            yp=-x_jd;
        end
        
        
        % 先找异常值，若有，看是否附近有超过阈值的Spike
        [~,L2] = findpeaks(yp,'minpeakheight',20);%x的rms=1
        [~,L] = findpeaks(yp,'minpeakheight',high,'MINPEAKDISTANCE',peak_interval);
        S(1,L)=1;
        if ~isempty(L2)
            % 有超过阈值的Spike，用临近点代替
            [~,L3] = findpeaks(yp,'minpeakheight',high);
            for k=1:length(L2)
                temp=abs(L3-L2(k));
                [num,loc]=min(temp);
                if num(1)<peak_interval/2 % 阈值<peak_interval
                    L(L==L2(k))=L3(loc(1));
                end
            end
            S(1,L)=1;
        end
        
        %         if high>0
        %             [~,L] = findpeaks(x,'minpeakheight',high,'MINPEAKDISTANCE',peak_interval);
        %             [~,L2] = findpeaks(x,'minpeakheight',param.out,'MINPEAKDISTANCE',peak_interval);%x的rms=1
        %             L3=setdiff(L,L2);%返回在L1中，不在L2中的元素
        %             S(1,L3)=1;
        %         elseif high<0
        %             [~,L] = median(-x,'minpeakheight',-high,'MINPEAKDISTANCE',peak_interval);
        %             [~,L2] = findpeaks(-x,'minpeakheight',param.out,'MINPEAKDISTANCE',peak_interval);%x的rms=1
        %             L3=setdiff(L,L2);%返回在L1中，不在L2中的元素
        %             S(1,L3)=1;
        % %             S(1,L)=1;
        %         else S=[];
        %         end
    end

end
