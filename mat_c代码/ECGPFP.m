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
    param.peakinterval=0.33*param.fs;
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
% x=single(x);%signal四个字节，节省空间。
if isempty(Spike) %判断输入是否为空
    residue=x;
else
    residue=x;
%     residue=bopi(x,Spike,param.wavelength,0);
end

% ---- parameter ----
count=1;

SpikeC=[];
W=[];
huifu_ECG=[];
huifu=[];

% yc=zeros(param.ficanum*param.outnum,size(x,2));%剥皮前FastICA运行数*输出数，信号采样数
source=[];
%tongdao=randperm(size(x,1));
%xw=yanchibaihua3(x(tongdao(1:min(size(x,1),64)),:),param.dc,0.0);
% xw=yanchibaihua3(x,param.dc,0.0);

% ---- starting  ----
while count<=param.peelnum
    % ---- delayed whiten  ----将SVD后的residue信号进行白化
%     ran=randi(param.df)
    [z,~,Cond]=yanchibaihua3(residue,param.df ,0);
    if Cond>1e15
        disp(['warn:the Cond is too large!',num2str(Cond)]);
        return
    end
    
    % ----deflation FastICA ----
    for i=1:param.ficanum
        [Se,Ae,cerr]=coshFpDeIca(z(1:floor((0.75+0.25*(rand>0.75))*size(z,1)),:),param.outnum);
        source=[source;Se];
%         yc(((i-1)*param.outnum+1):i*param.outnum,:)=Se;
        
        %         if param.debug
        %         multiPlot('FastICA',Se);
        %         end
        
    end
    % delete same yc
    yc=dis_spike(source,0.5);
    
    
    % ----threshold selecting ,clustering----
    F=[];
    for i=1:size(yc,1)
        [F1,F2,yc(i,:),th]=getspike(yc(i,:),param.fs,param.peakinterval,param.debug,ann);
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
    
    % delete same spike
%     jj=1;
%     if ~isempty(Spike)
%         while jj<=size(refer,1)
%         if cmp_spike(Spike,refer(jj,:),0.2)
%             refer(jj,:)=[];
%         else
%             jj=jj+1;
%         end
%         end
%     end
    
    % ----constraint fastICA ----
    xw=yanchibaihua3(x,param.dc,0.0);
    i=1;
    while i<=size(refer,1)
%         yc(i,:)=fcica3(xw,refer(i,:));
        figure
        plot(refer(i,:))
        hold on
        plot(ann,0.9,'r+')
        [yt, ~] = cfICA(xw, refer(i,:));
        ycmax = max(yt);
        figure;set(gcf,'Position',get(0,'ScreenSize'));
        hold on;
        plot(yt);title([num2str(i),'/',num2str(size(yt,1)),'Please input the threshold(or press "enter" to auto setting or input "0" to delete it):']);
        plot(ann,1,'r+');
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
            TH=spikethresh(yt);
            S=searchspike_in(yt,TH,param.peakinterval);
%             yt=fcica3(xw,S);
            [yt, ~] = cfICA(xw, S);
        else
            S=searchspike_in(yt,de,param.peakinterval);
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
                close all
                %开始下一轮FastICA
                break
                
                % ————判断与之前的spike是否重复————
            elseif ~cmp_spike(Spike,S,0.2)                
                %————加入序列————
                Spike=[Spike;S];
                
                %————STA 波形估计MECG,FECG————
                [residue,huifu,huifu_ECG,W]=bopi(x,Spike,param.wavelength,1);
                
                
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
                
                end
                return
            else
                close all
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
            flag=1;
        else
            flag=0;
        end
        
    end



end
