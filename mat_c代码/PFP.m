function [SpikeC,W,residue,huifu_ECG,huifu,peelcount,numc]=PFP(x,Spike,param)
if nargin==0
    disp('$$$ Documentation of Auto Progressive FastICA Peel-off (APFP)  $$$');
    disp('Author:Maoqi Chen ');
    disp('Email:hiei@mail.ustc.edu.cn ');
    disp('Update:2018.4.16');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp('# Form of the Function: Spike_Trains=PFP(EMG,Identified_Spike_Trains,param)');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp('# Input Arguments:');
    disp('EMG : the EMG signals (an n*T matrix,n channels, T samples).');
    disp('Identified_Spike_Trains : the identified spike trains (an p*T matrix, p identified MU numbers, T samples), the default value is [].');
    disp('param : a struct variable contains the important parameters of the algorithm. ');
    disp('        (%PS:The names of parameters in "param" should be the same as this documentation.%)')
    disp('        peelnum... the number of "peel-off" step, the default value is 3.');
    disp('        ficanum... the number of runs of FastICA before each "peel off" step, the default value is 5.');
    disp('        outnum...  the output number of each FastICA, the default value is 5.');
    disp('        fs... the sampling rate of EMG signal, the default value is 2kHz.');
    disp('        wavelength... the wavelength of the MUAPs, the default value is fs/20 samples. ');
    disp('        peakinterval... the minimum interval of successive spikes of one MU, the default value is fs/50.');
    
    disp('        df... the maximun delay factor of FastICA, the default value is 5 samples.');
    disp('        dc... the delay factor of constrained FastICA, the default value is 20 samples.');
    disp('        CORR... the CORRelation constraint between the output and the reference spike train in constrained fastica.');
    disp('                only the CORRelation coefficient is greater than CORR, we consider the output is a reliable MU spike train.');
    disp('                the default value of CORR is 0.5.');
    disp('        covs... the coefficient of variation (CoV) of the inter-spike intervals.');
    disp('                only the coefficient of variation is less than covs, we consider the output is a reliable MU spike train.');
    disp('                the default value of covs is 0.4.');
    disp('        covp... the coefficient of variation (CoV) of amplitudes of the output of constrained FastICA.');
    disp('                only the CORRelation coefficient is less than covp, we consider the output is a reliable MU spike train.');
    disp('                the default value of covs  is 0.35.');
    disp('        fig... whether to display the process of decomposition.');
    disp('                the default value of fig is 1 (to show the progress).');
    disp('        auto... automatic decomposition switch.');
    disp('                the default value of auto is 1 (automatic decomposition).');
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp('# some examples');
    disp('1. Spike_Trains=APFP(EMG);(use the default settings)');
    disp('2. Spike_Trains=APFP(EMG,Identified_Spike_Trains,param);(use the default settings in the case that we have got spike trains of some MUs)');
    disp('3. param.peelnum=10;Spike_Trains=APFP(EMG,[],param);(set the number of runs of "peel off" steps as 10)');
    disp('4. param.fig=0;Spike_Trains=APFP(EMG,[],param);(do not to display the process of decomposition)');
    disp('5. param.auto=0;Spike_Trains=APFP(EMG,[],param);(an interactive approach to set thresholds and select motor units)');
    return;
end %disp用来展示变量的内容，可以是字符串，元胞，矩阵，结构体。disp功能类似于c语言中的print;java语言中的System.out.println（)，
%可以输出几乎任何类型的变量。


% close all;
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
    param.peakinterval=60/230*param.fs;
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
x=x-mean(x,2)*ones(1,size(x,2));%x中每个元素减去均值
x=single(x);%signal四个字节，节省空间。
if isempty(Spike) %判断输入是否为空
    residue=x;
else
    residue=bopi(x,Spike,param.wavelength,0);
end
SpikeC=[];W=[];
count=1;
yc=zeros(param.ficanum*param.outnum,size(x,2));%剥皮前FastICA运行数*输出数，信号采样数
%tongdao=randperm(size(x,1));
%xw=yanchibaihua3(x(tongdao(1:min(size(x,1),64)),:),param.dc,0.0);
xw=yanchibaihua3(x,param.dc,0.0);
check=[];
CHECK=[];cache=[];peelcount=[];

if param.auto==1%auto
    while count<=param.peelnum
        if param.ini==1     %初始化参数
            a1=randi(size(residue,1),1,1);%生成一个1-size(residue,1)的1*1的矩阵
            a2=neo(residue(a1,:),1);%非线性能量算子
            a2=a2/norm(a2);          %norm(a)与norm(a,2)相同，向量所有值的平方和再开根号，单位化
            a3=spikethresh(a2);
            
            %~对应峰值，Lini对应峰值位数
            [~,Lini] = findpeaks(a2,'minpeakheight',a3,'MINPEAKDISTANCE',param.peakinterval);
        end
        
        [z,~,Cond]=yanchibaihua3(residue, randi(param.df)-1,0);%z为正交矩阵
        if Cond>1e10
            disp(['warn:the Cond is too large!',num2str(Cond)]);
            return
        end
        %tongdao=randperm(size(x,1));
        %xr=yanchibaihua3(residual(tongdao(1:min(size(x,1),64)),:),param.df,max(0.3,1-count/param.peelnum));
        %xr=yanchibaihua3(residual(tongdao(1:min(size(x,1),64)),:),param.df,0);
        %xr=yanchibaihua3(residual,param.df,0);
        for i=1:param.ficanum
            if param.ini==1
                a1=randi(size(residue,1),1,1);
                a2=neo(residue(a1,:),1); %非线性能量算子，将输出变为正的发放序列
                a2=a2/norm(a2);
                a3=spikethresh(a2); %得到自动迭代的阈值
                [~,Lini] = findpeaks(a2,'minpeakheight',a3,'MINPEAKDISTANCE',param.peakinterval);
                %a4是长度等于Lini的整数随机组合而成的向量
                a4=randperm(length(Lini));  %p = randperm(n) 返回行向量，其中包含从 1 到 n 没有重复元素的整数随机排列。
                Wini=Lini(a4(1:param.outnum));
                yc(((i-1)*param.outnum+1):i*param.outnum,:)=fastic_parallel3(z(1:floor((0.5+0.5*(rand>0.5))*size(z,1)),:),param.outnum,Wini);
                %运行fastica算法
            else
                yc(((i-1)*param.outnum+1):i*param.outnum,:)=fastic_parallel3(z(1:floor((0.5+0.5*(rand>0.5))*size(z,1)),:),param.outnum);
            end
            
            %yy(((i-1)*outnum+1):i*outnum,:)=fastic_parallel3(residual,randi(10)-1,outnum,0.8);
            
            %yc(((i-1)*param.outnum+1):i*param.outnum,:)=fastic_parallel3(residual,param.df,param.outnum,0.5*(rand>0.5));
        end
        %     [~,yc]=dis_spike(yy,peak_interval,1);
        %     if isempty(yc)
        %         count=count+1
        %         continue;
        %     end
        if ~isreal(yc)%isreal判断是不是实数
            disp('warn:output is not a real signal');
            return;
        end
        F=[];
        
        
        %去除重复的发放序列
        yc=dis_spike(yc,coef);%y_ini=initialspike(yc,param.peakinterval);
        
        for i=1:size(yc,1)
            [F1,F2,yc(i,:)]=getspike(yc(i,:),1,param.fs);temp=max(abs(yc(i,:)));
            if param.fig==1
                figure;set(gcf,'Position',get(0,'ScreenSize'));
                plot(yc(i,:));hold on;
                for jjj=1:size(F2,2)
                    scatter(F2{jjj},(-temp-2*jjj)*ones(size(F2{jjj})),'+');hold on;
                end
                
                pause(2);
                close
            end
            F=[F,F1];
        end
        if ~isempty(F)
            % 生成FastICA的0-1序列
            yc=zeros(size(F,2),size(xw,2));
            for i=1:size(F,2)
                yc(i,F{i})=1;
            end
        end
        
        
        
        i=1;check=[];
        while i<=size(yc)
            [S,~,CORR,cov_s,cov_p,much]=AutoCICAtest(xw,yc(i,:),param.peakinterval,param.fs,param.fig);
            nnaa=isnan(CORR*cov_s*cov_p);
            judgement1=(nnaa>0)|(CORR<param.CORR)|(cov_p>param.covp);
            if judgement1>0
                S=[];
                disp([CORR cov_p]);
            else
                judgement2=(cov_s>param.covs)|(much==0);
                if judgement2>0
                    check=[check;S];
                    S=[];
                end
            end
            Spike=[Spike;S];
            i=i+1;
        end
        Spike=dis_spike(Spike);
        check=dis_spike(check);
        
        %     figure;set(gcf,'Position',get(0,'ScreenSize'));
        %     for i=1:size(S,1)
        %         plot(S(i,:)-3*i);hold on;
        %     end;
        %     pause(2);close
        
        Cnum1=size(cache,1);
        [CHECK,cache]=Checkandcache(check,CHECK,cache);
        CHECK=dis_spike(CHECK);
        cache=dis_spike(cache);
        Cnum2=size(cache,1);
        SpikeC=[Spike;cache];SpikeC=dis_spike(SpikeC);peelcount(count)=size(SpikeC,1);
        if isempty(Spike)&&(Cnum2==Cnum1)
            count=count+1;
            continue;
        end
        if count==param.peelnum
            [residue,~,W]=bopi(x,SpikeC,param.wavelength,1);%取1会画图
        elseif (param.fig==1)
            residue=bopi(x,SpikeC,param.wavelength,1);close
        else residue=bopi(x,SpikeC,param.wavelength,0);
        end
        count=count+1;
    end
    %[residue,~,W]=bopi(x,SpikeC,param.wavelength,1);
    numc=size(Spike,1);
    
    
    
    %以下是手动分解步骤
else
    
    while count<=param.peelnum
        [z,~,Cond]=yanchibaihua3(residue, randi(param.df),0);
        if Cond>1e15
            disp(['warn:the Cond is too large!',num2str(Cond)]);
            return
        end
        for i=1:param.ficanum
            if param.ini==1
                a1=randi(size(residue,1),1,1);
                a2=neo(residue(a1,:),1);
                a2=a2/max(abs(a2));
                a3=spikethresh(a2);
                [~,Lini] = findpeaks(a2,'minpeakheight',a3,'MINPEAKDISTANCE',param.peakinterval);
                a4=randperm(length(Lini));
                Wini=Lini(a4(1:param.outnum));
                figure;subplot(2,1,1);plot(residue(a1,:));
                subplot(2,1,2);plot(a2);hold on;scatter(Wini,ones(param.outnum,1));pause(5);close
                
                
                yc(((i-1)*param.outnum+1):i*param.outnum,:)=fastic_parallel3(z(1:floor((0.5+0.5*(rand>0.5))*size(z,1)),:),param.outnum,Wini);
            else
                %yy(((i-1)*outnum+1):i*outnum,:)=fastic_parallel3(residual,randi(10)-1,outnum,0.8);
                % 将延迟白化信号的前面所有信号或一半信号放入FastICA中分解
                yc(((i-1)*param.outnum+1):i*param.outnum,:)=fastic_parallel3(z(1:floor((0.5+0.5*(rand>0.5))*size(z,1)),:),param.outnum);
                %yc(((i-1)*param.outnum+1):i*param.outnum,:)=fastic_parallel3(residual,param.df,param.outnum,0.5*(rand>0.5));
            end
        end
        multiPlot('FastICA',yc);
        
        if ~isreal(yc)
            disp('warn:output is not a real signal');
            return;
        end
        F=[];
        btw=[];
        yc=dis_spike(yc,coef);%y_ini=initialspike(yc,param.peakinterval);
        sname=strcat( num2str(param.user),'pfp_yc.mat');
        load(sname);
        for i=1:size(yc,1)
%             需要聚类，因为存在噪声时，Spike可能包含噪声Spike和心电Spike
            [F1,F2,yc(i,:),th]=getspike(yc(i,:),1,param.fs);temp=max(abs(yc(i,:)));
            if param.fig==1
                figure;set(gcf,'Position',get(0,'ScreenSize'));
                plot(yc(i,:));hold on;
                plot(th*ones(1,size(yc,2)),'r-');
%                 yline(th,'r--');hold on;
                for jjj=1:size(F2,2)
                    scatter(F2{jjj},(temp/2-1*jjj)*ones(size(F2{jjj})),'+');hold on;
                end
%                 for jjj=1:size(F1,2)
%                     scatter(F1{jjj},(-temp-2*jjj)*ones(size(F1{jjj})),'+');hold on;
%                 end
%                                 pause(2);
%                                 close
            end
            
            F=[F,F1];
        end
        if ~isempty(F)
            yc=zeros(size(F,2),size(xw,2));
            for i=1:size(F,2)
                if ~isempty(F{i})
                    yc(i,F{i})=1;
                end
            end
        end
        
        i=1;
        while i<=size(yc)
            yc(i,:)=fcica3(xw,yc(i,:));
            %             if skewness(yc(i,:)')<0
            %                 yc(i,:)=-yc(i,:);
            %             end
            ycmax = max(yc(i,:));
            figure;set(gcf,'Position',get(0,'ScreenSize'));
            plot(yc(i,:));title([num2str(i),'/',num2str(size(yc,1)),'Please input the threshold(or press "enter" to auto setting or input "0" to delete it):']);
            
            %此循环确保de为数字，之前没有循环有可能会误触导致程序中断
            %此外如果输入字母i，j不会认为是字母，会认为是复数的虚数单位，故要首先排除i，j
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
            if isempty(de)
                TH=spikethresh(yc(i,:));
                S=searchspike_in(yc(i,:),TH,param.peakinterval);
                yt=fcica3(xw,S);
            else
                S=searchspike_in(yc(i,:),de,param.peakinterval);
                yt=fcica3(xw,S);
            end
            ytmax = max(yt);
            figure;set(gcf,'Position',get(0,'ScreenSize'));
            subplot(211);plot(yt);title('Determine the threshold (or "enter" to accept or "0" to discard):');subplot(212);plot(S);
            
            %第一次确定阈值，接受或者重新输入阈值，此时enter变为接受，只能手动输入
            %此外如果输入字母i，j不会认为是字母，会认为是复数的虚数单位，故要首先排除i，j
            while(1)
                de=input('Determine the threshold (or "enter" to accept or "0" to discard):','s');
                if strcmp(de ,'i')||strcmp(de ,'j')
                    continue;
                end
                if (isempty(de))
                    break;
                else
                    de = str2num(de);
                    if (~isempty(de))&& (de<ytmax)
                        break;
                    end
                end
            end
            
            close
            while (~isempty(de) & de~=0)
                S=searchspike_in(yt,de,param.peakinterval);
                yt=fcica3(xw,S);
                ytmax = max(yt);
                figure;set(gcf,'Position',get(0,'ScreenSize'));
                subplot(211);plot(yt);title('Determine the threshold (or "enter" to accept or "0" to discard):');subplot(212);plot(S);
                
                %重新设置阈值直到删除/接受
                %此外如果输入字母i，j不会认为是字母，会认为是复数的虚数单位，故要首先排除i，j
                while(1)
                    de=input('Determine the threshold (or "enter" to accept or "0" to discard):','s');
                    if strcmp(de ,'i')||strcmp(de ,'j')
                        continue;
                    end
                    if (isempty(de))
                        break;
                    else
                        de = str2num(de);
                        if (~isempty(de))&& (de<ytmax)
                            break;
                        end
                    end
                end
                
                close
            end
            if de==0
                S=[];
            end           
            Spike=[Spike;S];i=i+1;
        end
         % 比较是否出现重复序列
        Spike=dis_spike(Spike,0.2);peelcount(count)=size(Spike,1);
        if isempty(Spike)
            count=count+1;
            continue;
        end
        
        % 排序
        temp=[];
        for j=1:size(Spike,1)
            temp(j)=length(find(Spike(j,:)==1));
        end
        if size(Spike,1)>=2
            if temp(1)>temp(2)
                temp_S=Spike(1,:);
                Spike(1,:)=Spike(2,:);
                Spike(2,:)=temp_S;
            end
        end
        
        [residue,huifu,huifu_ECG,W]=bopi(x,Spike,param.wavelength,1);title('Press "enter" to continue or "0" to finish):')
        re=input('Press "enter" to continue or "0" to finish):');
        if ~isempty(re)
            %             figure;
            %             for i=1:size(Spike,1)
            %                 plot(Spike(i,:)-3*i);hold on;
            %             end
            title('')
            SpikeC=Spike;numc=size(SpikeC,1);
            return;
        end
        if count==param.peelnum
            title('');
        else
            close
        end
        count=count+1;
    end
    %figure;set(gcf,'Position',get(0,'ScreenSize'));
    SpikeC=Spike;numc=size(SpikeC,1);
    %[residue,~,W]=bopi(x,Spike,param.wavelength,1);%2018/2/24注释掉的
    %     figure;
    %     for i=1:size(Spike,1)
    %         plot(Spike(i,:)-3*i);hold on;
    %     end
end

% figure;set(gcf,'Position',get(0,'ScreenSize'));
% for i=1:size(SpikeC,1)
%     plot(SpikeC(i,:)-3*i);hold on;
% end


%剥皮
    function [residual,huifu,huifu_ECG,W]=bopi(y,S,wavelength,p)%pwei1zehuatu
        a=randi(size(y,1),1,1);
        if isempty(S)
            residual=y;
            W=[];
            huifu=y;
            return
        end
        [~,xt,~,huifu_ECG,C]=boxing(y(a,:),S,wavelength,p);
        huifu=y*C'*xt;
        residual=y-huifu;
        W=C*y';
        
        function [W,xt,huifu,huifu_ECG,C]=boxing(x,S,wavelength,p)%x行向量。spike矩阵
            %% STA波形估计 LS
            % 动态确定波形长度
            [n,T]=size(S);%n通道数，t样本点数？            
            yanchi=fix(0.45*wavelength);
            xt=zeros(n*wavelength,T);xt=single(xt);S=single(S);%xt是一个扩展矩阵
            
            for j=1:n
                %托普利兹矩阵
                xt((j-1)*wavelength+1:j*wavelength,:)=toeplitz([S(j,(yanchi+1):-1:1),zeros(1,(wavelength-yanchi-1))],[S(j,(yanchi+1):end),zeros(1,yanchi)]);
            end
            C=(xt*xt')\xt;
            W=C*x';
            huifu=W'*xt;
            
%             xls = pinv(x_in'*x_in)*x_in'*y_in;
%             % Generate LAD coeff 
%             % Setup reformulated LP variables
%             A = xt';
%             b = x';
%             
%             % LAD-TO-LP REFORMULATION with suppression trick
%             len_x = size(A,2);
%             c = [zeros(len_x,1);ones(size(b))];
%             F = [A -eye(size(b,1)); -A -eye(size(b,1))];
%             g = [b; -b];
%             %F = [A -eye(size(b,1)); -A -eye(size(b,1));zeros(size(A)) -eye(size(b,1))];
%             %g = [b ; -b; zeros(length(b),1)];
%             
%             % Run the LP solver
%             z = linprog(c,F,g);
%             xlad = z(1,:);
            
            %在此处进行ECG绘制
            huifu_ECG=[];
            for j=1:n
                huifu_ECG=[huifu_ECG;W((j-1)*wavelength+1:j*wavelength,:)'*xt((j-1)*wavelength+1:j*wavelength,:)];                
            end
                        
            if (p==1)
                figure;set(gcf,'Position',get(0,'ScreenSize'));
                maxi=1.5*max(max(abs(x)),max(abs(huifu)));
                subplot(2,4,2:4);plot(x);hold on;
                for j=1:n
                    plot(huifu_ECG(j,:)-j*maxi);hold on; 
%                     plot(huifu_ECG(:,:,j)-j*maxi);hold on;                
                end
                plot(huifu-(j+1)*maxi);hold on;
                plot(x-huifu-(j+2)*maxi);hold on;
                subplot(2,4,6:8);
                for ii=1:size(S,1)
                    plot(S(ii,:)-2*ii);hold on;
                end
                hold on;
                subplot(n,4,1:4:4*n-3);jg=max(abs(W));
                for j=1:n
%                     plot(Ws(:,:,j)-1*jg*j,'linewidth',3);hold on;
                    plot(W(1+wavelength*(j-1):wavelength*j)-2*jg*j,'linewidth',3);hold on;
                end
                pause(2);
            end
        end
        
    end

%白化去相关
% 这里白化的意思是：
% 对原始信号X进行预处理 预处理包括去均值和白化(Whitening)，是通过对观测数据向量进行线性变换，使其均值为零，方差为1，去除各观测之间的相关性。
%
% 白化过程 简单而言就是 将信号或者噪声的协方差矩阵的对角化处理，不同的信号或信号形式其协方差矩阵一般而言是不一样的，那么我们从矩阵论的知识可以得到，两者要求的正交矩阵也是不一样的。
%R是the delay factor of constrained FastICA, the default value is 20 samples
%x为原始信号 Nsignal通道数
    function [y,Q,Cond] = yanchibaihua3(x, R,po)
        % function xt= yanchibaihua(x, R)
        [Nsignal, Nsample] = size(x);
        xt = zeros((2*R+1)*Nsignal,Nsample);
        for k= 1:Nsignal
            %xt(R*(i-1)+1:R*i,:)=toeplitz([x(i,1);zeros(R-1,1)],x(i,:)); toeplitz(x)用向量x生成一个对称的托普利兹矩阵
            xt((2*R+1)*(k-1)+1:(2*R+1)*k,:)=toeplitz([x(k,(R+1):-1:1),zeros(1,R)],[x(k,(R+1):end),zeros(1,R)]);%xt为扩展矩阵
        end
        xt=xt-mean(xt,2)*ones(1,Nsample);%将xt中心化 ,求每行的均值
        %  mxt = xt*xt'/Nsample;
        x_cov=cov(xt');                    % cov为求协方差的函数
        %当Cond很大时，可认为x_cov为奇异矩阵，即行列式为0.
        Cond=cond(x_cov);%矩阵A的条件数等于A的范数与A的逆的范数的乘积，即cond(A)=‖A‖·‖A^(-1)‖，对应矩阵的3种范数，相应地可以定义3种条件数。 函数 cond(A,1)、cond(A)或cond(A inf) 是判断矩阵病态与否的一种度量，条件数越大矩阵越病态。
        [E,D]=eig(x_cov);  % 对信号矩阵的协方差函数进行特征值分解 E是特征向量 D是特征值(只有对角线有值)
        % [E,D]=eig(mxt);
        lamda=diag(D); %diag是(提取对角元素)
        %将特征根及其特征向量按降序的方式重新排列
        num=length(lamda);trun=fix(num*(1-po));
        [w,j]=sort(lamda,'descend');E=E(:,j(1:trun)');D=diag(w(1:trun));%D排序好了，E每行顺序调换后和D顺序一致
        
        
        Q=sqrt(D)\(E)';                        % Q为白化矩阵
        y=Q*xt;                      %y为正交矩阵    MixedS_white为白化后的信号矩阵
        % IsI=cov(y');                   %近似单位阵
    end
%白化


%FASTICA
    function [y_out,w]=fastic_parallel3(z,num,wini)%这里的w正好和fastic_parallel2中的w互为转置，num wini
        %[z,~]=yanchibaihua3(x, R,po);%z是正交矩阵
        [m,T]=size(z);
        if nargin<3  %nargin为"number of input arguments"的缩写。 在matlab中定义一个函数时， 在函数体内部， nargin是用来判断输入变量个数的函数。
            w=randn(num,m);
            w=sqrtm(w*w')\w;%randn：产生均值为0，方差σ^2 = 1，标准差σ = 1的正态分布的随机数或矩阵的函数。sqrtm求矩阵的平方根
        else
            w=z(:,wini)';w=sqrtm(w*w')\w; %w初始权矢量（任意）
        end
        w1=randn(num,m);w1=sqrtm(w1*w1')\w1;%左除式A\B，则相当于inv（A）*B
        
        c1=1;%c1迭代次数
        while (1-min(abs(diag(w*w1')))>1e-5)&&(c1<1000)%0.0000000000001 %判断是否收敛
            w1=w;
            w=g(w1*z)*z'-diag(gd(w1*z)*ones(T,1))*w1;
            w=sqrtm(w*w')\w;
            c1=c1+1;
            %     [E,D]=eig(w*w');
            %     w=E/(sqrt(D))*E'*w;
        end
        y_out=w*z; % y_out输出的独立分量
        if c1>=1000
            disp('warn:FASTICA does not converged whin 1000 iterations.');
        end
        
        function G=g(x) %G为非线性函数,g(x)为导数,tanh(x)为log(cosh(x))的导数
            G=tanh(x);
            
        end
        function G=gd(x) %gd(x)为二次导数
            G=1-(tanh(x)).^2;
            
        end
    end
%FASTICA

%区分spike，根据相关性 删除重复的发放序列
    function [S]=dis_spike(y,p,fs)
        if nargin==1
            p=0.4;
        end
        if nargin==2
            fs=1000;
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
        multiPlot('shaixuanguohoudey',S);
    end




%CICA   带有时域约束的fastICA
    function [S,yt,xg,coef_s,coef_w,much]=AutoCICAtest(xw,y,peak_interval,fs,pic)
        %              yp=y;yp(yp<0)=0;yn=-y;yn(yn<0)=0;
        %             th1=spikethresh(yp);
        %             th2=spikethresh(yn);
        %             if th1>th2
        %             S=searchspike_in(yp,th1,peak_interval);
        %             else
        %             S=searchspike_in(yn,th2,peak_interval);
        %             end
        %         if skewness(y')<0   %skewness返回向量的偏斜度，看众数是在算数平均数的左侧还是右侧
        %             y=-y;
        %         end
        th2=spikethresh(y);
        S=searchspike_in(y,th2,peak_interval);
        %%这个是debug用以后删掉
        %         if isempty(S)
        %             figure;plot(y);
        %             title('bug');
        %             pause
        %         end
        %%%%
        [yt,xg]=fcica3(xw,S);
        th2=spikethresh(yt);th1=0;c=1;
        if pic==1
            figure;set(gcf,'Position',get(0,'ScreenSize'));
            while abs(th1-th2)>1e-3 && c<20
                th1=th2;
                S=searchspike_in(yt,th1,peak_interval);
                [yt,xg]=fcica3(xw,S);
                %th2=spikethresh(yt);
                [th2,lo]=spikethresh(yt);
                plot(yt);hold on;plot(th2*ones(size(yt)));hold on;plot(lo*ones(size(yt)));pause(1);hold off;
                c=c+1;
            end
            close
        else
            while abs(th1-th2)>1e-3 && c<20
                th1=th2;
                S=searchspike_in(yt,th1,peak_interval);
                [yt,xg]=fcica3(xw,S);th2=spikethresh(yt);
                c=c+1;
            end
        end
        % 波形变异系数
        a=diff(find(S==1));
        much=(length(a(a<fs/2))>=sum(a(a<fs/2))/fs*3);%发放个数的限制
        a((a>fs/4))=[];
        b=yt(S==1);
        %         if numel(b)>4
        %         lillietest(b)
        %         end
        coef_s=std(a)/mean(a);
        coef_w=std(b)/abs(mean(b));
        %         xg
        if c==30
            disp('!!Warn:the automatic threshold does not converge');
        end
    end


    function [y,xg,y0,miu,sigma,ksi]=fcica3(xw,S)%约束q=(ksi^4-E{(yr)^4})/4
        %x=yanchibaihua2(x,10);
        %x=yanchibaihua3(x,15,0);
        T=size(xw,2);
        %r=baihua(S);
        r=S/sqrt(mean(S.^2));
        sigma=1e1;ksi=1; thr=1e-5;
        flag=1;c=40;
        figure
        while flag
            sigma=1e1;ksi=1;
            %w1=rand(size(xw,1),1);w1=w1/norm(w1);
            w1=zeros(size(xw,1),1);
            w2=mean(xw(:,S>0),2);w2=w2/norm(w2);
            y0=w2'*xw;
            % gaosi=randn(1,size(xw,2));
            % Ggs=mean(G(gaosi));
            miu=0;k=1;
            aa=abs(abs(w1'*w2)-1);
            
            while aa>thr
                sigma=sigma*1.02;
                w1=w2;
                y=w1'*xw;
                w2=1/T*xw*g(y)'+miu/T*xw*r';
                w2=w2/norm(w2);
                aa=abs(abs(w1'*w2)-1);
                xg=mean((w2'*xw.*r));
                miu=max([0,miu+sigma*(ksi-xg)]);
                if ksi>=xg
                    ksi=ksi*0.97;
                end
                k=k+1;
                if k>c
                    thr=thr*10;
                    break;
                end
            end
            flag=0;
            if k>c
                c=c+10;
                flag=1;
            end
            win=5;
             yy=w2'*xw;
            plot(yy(2*win+1:end-2*win));
        end
        w=w2;
        y=w'*xw;
        % figure;subplot(212);plot(y);subplot(211);plot(r);
        %         function G=G(x)
        %             %G=log(cosh(x));
        %             G=x.*exp(-x.^2/2);
        %         end
        
        function G=g(x)
            % G=tanh(x);
            G=x.*exp(-x.^2/2);
        end
        %         function G=gd(x)
        %             %G=1-(tanh(x)).^2;
        %             G=(1-x.^2).*exp(-x.^2/2);
        %         end
    end


    function Neo=neo(data,n)   %计算非线性能量算子的函数
        [m,T]=size(data);
        data=data-mean(data,2)*ones(1,T);
        Neo=zeros(m,T);
        for j=n+1:T-n
            Neo(:,j)=(data(:,j)).^2-data(:,j-n).*data(:,j+n);
        end
        % figure;
        % for i=1:m
        %     subplot(m,1,i);plot(Neo(i,:));
        % end
    end


%计算自动迭代阈值
    function [th,th0]=spikethresh(x,p)
        if nargin<2
            p=0.5;
        end
        if all((x==1)==x)
            th=0.5;
            return;
        end
        H=median(abs(x))/0.6754;%median求中位数
        if (H>=max(abs(x)))||(H==0)
            H=0.001;
        end
        %x(x<=H)=H;
        x=x(x>=H);x(x>0.95*max(x))=0.95*max(x);
        %a=sort(x,'descend');
        TT=size(x,2);
        seg=201;
        tth=max(x)/seg:max(x)/seg:(seg-1)*max(x)/seg;g=zeros(1,seg-1);
        for q=1:(seg-1)
            p0=length(x(x>=tth(q)))/TT;
            g(q)=p0*(1-p0)*(mean(x(x>tth(q)))-mean(x(x<=tth(q)))).^2;
        end
        lo=find(g==max(g), 1, 'first' );
        th=tth(lo);
        %th=mean(a(1:20))/2;
        th0=tth(lo);
        t1=10000;c=1;
        while abs(t1-th)>1e-3 && c<=100
            t1=th;
            th=mean(x(x>=t1))*p+mean(x(x<t1))*(1-p);
            %th=max(mean(x(x>=t1))*p+mean(x(x<t1))*(1-p),H);
            c=c+1;
        end
        %figure;set(gcf,'Position',get(0,'ScreenSize'));plot(x);hold on; plot(th*ones(size(x)));pause(2);close;
    end



    function [F,F2,x,th,Spike,Sc,L,C]=getspike(x,p,fs)%解决符号不一致的一些问题
%         interval=60/300*fs;%fix是靠0取整
        %         interval=fix(fs/100);%fix是靠0取整
        win=fix(fs/100);
        x_jd=x(2*win+1:end-2*win);
        yp=x_jd;yp(yp<0)=0;
        yn=-x_jd;yn(yn<0)=0;
        yp_fil=yp;
        yn_fil=yn;
        yp_fil(abs(yp)>20)=0;        
        yn_fil(abs(yn)>20)=0;
        
%         [~,L1] = findpeaks(yp,'minpeakheight',2.2,'MINPEAKDISTANCE',param.peakinterval);
%         [~,L2] = findpeaks(yp,'minpeakheight',18,'MINPEAKDISTANCE',param.peakinterval);%x的rms=1
%         L3=setdiff(L1,L2);%返回在L1中，不在L2中的元素
%         th1=sum(x_jd(L3));
%         
%         [~,L1i] = findpeaks(yn,'minpeakheight',2.2,'MINPEAKDISTANCE',param.peakinterval);
%         [~,L2i] = findpeaks(yn,'minpeakheight',18,'MINPEAKDISTANCE',param.peakinterval);%x的rms=1
%         L3i=setdiff(L1i,L2i);%返回在L1中，不在L2中的元素
%         th2=sum(x_jd(L3i))
%         
%         if abs(th1)<abs(th2)
%             x=-x;
%             L=L3i;
%         else
%             L=L3;
%         end
        
        th1=spikethresh(yp_fil)
        th2=spikethresh(yn_fil)
        if th1>th2
            % S=searchspike_in(yp,th1,peak_interval);
            [~,L1] = findpeaks(yp,'minpeakheight',th1,'MINPEAKDISTANCE',param.peakinterval);
%             [~,L1] = findpeaks(yp,'minpeakheight',th1,'MINPEAKDISTANCE',param.peakinterval);
            [~,L2] = findpeaks(yp,'minpeakheight',param.out,'MINPEAKDISTANCE',param.peakinterval);%x的rms=1
            L3=setdiff(L1,L2);%返回在L1中，不在L2中的元素
%             th1=mean(x_jd(L3));
            L=L3;
            th=th1;
        else
            % S=searchspike_in(yn,th2,peak_interval);
            [~,L1i] = findpeaks(yn,'minpeakheight',th2,'MINPEAKDISTANCE',param.peakinterval);
%             [~,L1i] = findpeaks(yn,'minpeakheight',th2,'MINPEAKDISTANCE',param.peakinterval);
            [~,L2i] = findpeaks(yn,'minpeakheight',param.out,'MINPEAKDISTANCE',param.peakinterval);%x的rms=1
            L3i=setdiff(L1i,L2i);%返回在L1中，不在L2中的元素
%             th2=mean(x_jd(L3i));
            L=L3i;
            x=-x;
            th=abs(th2);
        end
  
     
        L=L+2*win;
        if numel(L)<length(x(1,:))/param.fs/3  %返回L中元素的个数
            F=[];F2=[];Spike=[];Sc=[];L=[];C=[];
            return;
        end
        
         % 将x归一化
        x_n=x;
        x_n(x<0)=0;
        x_n = normalize(x_n)*10;
        
        Spike=zeros(length(L),2*win+1);
        for ii=1:length(L)
            Spike(ii,:)=x_n(L(ii)-win:L(ii)+win);
%             [~,pos]=max(x(L(ii)-win:L(ii)+win));
%             pos=pos+L(ii)-win-1;
%             Spike(ii,:)=x(pos-win:pos+win);
        end
          
        Sc(:,1)=sqrt(mean(Spike.^2,2));%每一行的平方根
        Sc(:,2)=sqrt(mean(diff(Spike,1,2).^2,2));%diff差分运算，1，2指每行后面的减去前面的
        %             figure;
        %             scatter(Sc(:,1),Sc(:,2));
        [C]=valley_seeking(Sc,20,0.3,3);
        if p==1
            figure;
            for ii=1:max(C)
                scatter(Sc(C==ii,1),Sc(C==ii,2)); hold on
            end
            scatter(Sc(C==-1,1),Sc(C==-1,2));hold on
%             pause(3);close
            legend('valley seeking');
        end
%         F{1}=L;
        F=[];
        F2=cell(1,max(C));
        index=1;
        for ii=1:max(C)
            F2{ii}=L(C==ii);
            if length(L(C==ii))>length(x(1,:))/param.fs/3
                F{index}=L(C==ii);
                index=index+1;
            end
        end
        
    end

    function y=normalize(x)
        s=max(abs(x));
        y=x/s;
    end

%找到x中非异常值的峰值的位置
    function S=searchspike_in(x,high,peak_interval,fs)
        if nargin==3
            fs=1000;
        end
        
        S=zeros(size(x));
        win=fix(fs/200);        
        x_jd=x(2*win+1:end-2*win);
        if high==0
            S=[];
            return
        elseif high>0
            yp=x_jd;
        else
            yp=-x_jd;
        end
        % 先找异常值，若有，看是否附近有超过阈值的Spike，代替
        [~,L2] = findpeaks(yp,'minpeakheight',50,'MINPEAKDISTANCE',peak_interval/2);%x的rms=1
        [~,L] = findpeaks(yp,'minpeakheight',high,'MINPEAKDISTANCE',peak_interval);
        if isempty(L2)
            S(1,L)=1;
        else
            [~,L3] = findpeaks(yp,'minpeakheight',high,'MINPEAKDISTANCE',peak_interval/2);
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



%寻谷聚类算法
    function [C,Cd,J]=valley_seeking(x,t1,t2,t3)
        %t1:局部临近，t3：局部凸，t2：聚类同质性
        D=squareform(pdist(x));d=size(D,1);%pdist(x)计算各行的距离，D为各点距离方阵，对角线元素为0
        [~,l]=sort(D,2);L=zeros(d,d);%~为D中每行按从小到大排列，l为相应位置的索引，L为NN
        for ii=1:d
            L(ii,l(ii,:))=0:d-1;
        end
        Ss=(L+L')/2;
        %J：NDD
        J=abs(L-L')./(Ss.^(1+1/d));
        [~,L_s]=sort(L,2);% x 的第 i 个邻居为L_S(i)
        J(isnan(J))=0;
        %观察J的分布

        D2=zeros(d,d);
        for ii=1:d
            for j=1:d
                if ii==j
                    D(ii,j)=0;
                    continue;
                end
                
                ka=L(ii,j)*sum(L(L_s(ii,1:L(ii,j)-1),ii));
                lb=L(j,ii)*sum(L(L_s(j,1:L(j,ii)-1),j));
                ki=L(ii,j)*sum(L(ii,L_s(ii,1:L(ii,j)-1)));
                lj=L(j,ii)*sum(L(j,L_s(j,1:L(j,ii)-1)));
                D2(ii,j)=(ka+lb)/(ki+lj);
%                 D2(ii,j)=(L(ii,j)*sum(L(L_s(ii,1:L(ii,j)),ii))+L(j,ii)*sum(L(L_s(j,1:L(j,ii)),j)))/(L(ii,j)*L(ii,j)*(L(ii,j)-1)+L(j,ii)*L(j,ii)*(L(j,ii)-1));
            end
        end
        D2(isnan(D2))=1;
        Cd=(Ss<=t1).*(J<=t2).*(D2<=(t3));
        %hist(J);pause;close
        Total=1:d;
        C=zeros(d,1);k=1;
        while sum(C==0)
            X=Total(C==0);
            A=X(1);
            B=find(sum(Cd(A,:),1));
            while (~isempty(setdiff(B,A)))
                A=union(A,B);
                B=find(sum(Cd(A,:),1));
            end
            if numel(A)<d/50
                C(A)=-1;
            else
                C(A)=k;
                k=k+1;
            end
        end
        
    end


%判断所得到的发放序列是否可靠
    function [CHECK,cache]=Checkandcache(check,CHECK,cache)
        if size(CHECK,1)==0
            CHECK=check;
            cache=[];
            return
        end
        La=zeros(size(check,1),1);
        for uu=1:size(check,1)
            Co=zeros(size(CHECK,1),1);
            for j=1:size(CHECK,1)
                Co(j)=max(xcorr(check(uu,:),CHECK(j,:),'coeff'));
            end
            if max(Co)>0.5
                La(uu)=1;
            end
        end
        cache=[cache;check(La==1,:)];
        CHECK=[CHECK;check(La==0,:)];
    end

%     function th=spikethresh(x,p)
%         if nargin<2
%             p=0.5;
%         end
%         if all((x==1)==x)
%             th=0.5;
%             return;
%         end
%         %th=median(abs(x));t1=0;c=1;
%         a=sort(x,'descend');H=median(abs(a))/0.6754;
%         if (H>=max(abs(x)))||(H==0);
%             H=0.001;
%         end
%         %x(x<=H)=H;
%         x=x(x>=H);
%         th=mean(a(10:20))/2;
%         if th<=H
%             th=H;return
%         end
%         t1=10000;c=1;
%         while abs(t1-th)>1e-3 && c<=100
%             t1=th;
%             th=mean(x(x>=t1))*p+mean(x(x<t1))*(1-p);
%             %th=max(mean(x(x>=t1))*p+mean(x(x<t1))*(1-p),H);
%             c=c+1;
%         end
%         %figure;set(gcf,'Position',get(0,'ScreenSize'));plot(x);hold on; plot(th*ones(size(x)));pause(2);close;
%     end
%     function [x_white]=baihua(x)
%         %去均值
%         x_mean=mean(x,2);
%         x_M=x_mean*ones(1,size(x,2));
%         x_new=x-x_M;
%
%         x_cov=cov(x_new');                    % cov为求协方差的函数
%         [E,D]=eig(x_cov);                      % 对信号矩阵的协方差函数进行特征值分解
%         Q=sqrt(D)\(E)';                        % Q为白化矩阵
%         x_white=Q*x_new;                      % MixedS_white为白化后的信号矩阵
%         % IsI=cov(x_white');                   %近似单位阵
%     end


%     function [y_out,w2]=fastica_single(z,F)
%         %[z,~]=yanchibaihua3(x, R,0);
%         [m,T]=size(z);
%         w1=randn(m,1);
%         w1=w1/norm(w1);
%         w2=mean(z(:,F),2);w2=w2/norm(w2);c=1;
%
%         while (abs(1-abs(w1'*w2))>1e-6)&&(c<2000)
%             w1=w2;
%             w2=(1/T*g(w1'*z)*z')'-mean(gd(w1'*z))*w1;w2=w2/norm(w2);c=c+1;
%         end
%         if c>=2000
%             abs(1-abs(w1'*w2))
%         end
%         %w=w2;
%         y_out=w2'*z;
%         %plot(y_out);pause;close;
%         %ww=w'*Q;
%         function G=g(x)
%             G=tanh(x);
%             %       G=x.*exp(-x.^2/2);
%             %       G=0.5*(1+sign(x).*(1-exp(-abs(x)/0.001)));%0.5*(1+sign(t).*(1-exp(-abs(t)/b))
%             %        G=x.*exp(-3.348*abs(x));
%             %        b=.5;G=x./(1+b*abs(x)).^2;
%         end
%         function G=gd(x)
%             G=1-(tanh(x)).^2;
%             %         G=(1-x.^2).*exp(-x.^2/2);
%             %        G=0.5/0.001*exp(-abs(x)/0.001);%0.5/b*exp(-abs(t)/b)
%             %          G=exp(-3.348*abs(x)).*(1-3.348*abs(x));
%             %       b=.5;G=(1-b*abs(x))./(1+b*abs(x)).^3;
%         end
%     end
% function label=dis_w(w)
%         A=1:size(w,2);label=[];
%         while ~isempty(w)
%             label=union(label,A(1));
%             W=abs(w'*w)>0.1;
%             a=find(W(:,1));
%             A(a)=[];
%             w(:,a)=[];
%         end
%     end
%     function [S,yt,CORR,coef_s,coef_w]=AutoCICAtest2(xw,S,peak_interval)
%         coef_s=zeros(1,size(S,1)); coef_w=zeros(1,size(S,1)); CORR=zeros(1,size(S,1));
%         yt=zeros(size(S,1),size(xw,2));
%         if size(S,1)>1
%             p=0.8;
%         else p=0.5;
%         end
%         for kk=1:size(S,1)
%             [yt(kk,:),xg]=fcica3(xw,S(kk,:));
%             th2=spikethresh(yt(kk,:),p);th1=0;c=1;
%             while abs(th1-th2)>1e-2 && c<10
%                 th1=th2;
%                 S(kk,:)=searchspike_in(yt(kk,:),th1,peak_interval);
%                 [yt(kk,:),xg]=fcica3(xw,S(kk,:));th2=spikethresh(yt(kk,:),p);c=c+1;
%             end
%             a=diff(find(S(kk,:)==1));a((a<40|a>500))=[];
%             b=yt(kk,S(kk,:)==1);
%             coef_s(kk)=std(a)/mean(a);
%             coef_w(kk)=std(b)/abs(mean(b));
%             CORR(kk)=xg;
%              subplot(211);plot(yt(kk,:));subplot(212);plot(S(kk,:));pause(10);close
%             if c==10
%                 disp('！Warn 自动设置的阈值未收敛');
%             end
%         end
% %             coef_s
% %             coef_w
% %             CORR
%     end


%     function [y_out,w]=fastica_parallel4(z,F)%这里的w正好和fastic_parallel2中的w互为转置
%         num=numel(F);
%         a=fix(size(z,1)*(1-0.5*(rand>0.5)));z=z(1:a,:);
%         [m,T]=size(z);
%
%         w1=randn(num,m);w1=sqrtm(w1*w1')\w1;
%         for cc=1:num
%             w(cc,:)=mean(z(:,F{cc}),2);
%         end
%         w=sqrtm(w*w')\w;
%         %         [~,weizhi]=sort(neo(z(randi(m),:),1),'descend');
%         %         A = randperm(50);B = A(1:num);
%         %         w=z(:,weizhi(B))';w=sqrtm(w*w')\w;
%         c1=1;
%         while (1-min(abs(diag(w*w1')))>1e-5)&&(c1<5000)%0.0000000000001
%             w1=w;
%             w=g(w1*z)*z'-diag(gd(w1*z)*ones(T,1))*w1;
%             w=sqrtm(w*w')\w;
%             c1=c1+1;
%             %     [E,D]=eig(w*w');
%             %     w=E/(sqrt(D))*E'*w;
%         end
%         y_out=w*z;
%         if c1>=5000
%             c1
%         end
%
%         function G=g(x)
%             G=tanh(x);
%
%         end
%         function G=gd(x)
%             G=1-(tanh(x)).^2;
%
%         end
%     end


%     function S=deleteyc(s)
%         S=[];
%         while ~ isempty(s)
%             S=[S;s(1,:)];
%             Co=zeros(size(s,1),1);
%             for j=1:size(s,1)
%                 Co(j)=max(xcorr(s(1,:),s(j,:),'coeff'));
%             end
%             s(Co>=0.4,:)=[];
%         end
%     end
% function [y_out,w]=fastica_1by1(z,F)%这里的w正好和fastic_parallel2中的w互为转置
%         num=numel(F);
%         %a=fix(size(z,1)*(1-0.1*(rand>0.5)));z=z(1:a,:);
%         %a=fix(size(z,1)*0.2);z=z(1:a,:);
%         %z=z(1:60,:);
%         [m,T]=size(z);
%         w=zeros(num,m);
%         %w1=randn(num,m);w1=sqrtm(w1*w1')\w1;
%         for cc=1:num
%             w(cc,:)=mean(z(:,F{cc}),2);w(cc,:)=w(cc,:)/norm(w(cc,:));
%         end
%         for k=1:num
%         c=1;w1=randn(1,m);
%         while (abs(1-abs(w(k,:)*w1'))>1e-6)&&(c<5000)
%     w1=w(k,:);
%     w(k,:)=1/T*g(w1*z)*z'-mean(gd(w1*z))*w1;w(k,:)=w(k,:)/norm(w(k,:));c=c+1;
%         end
%         if c>=5000
%             c
%         end
%         end
%         y_out=w*z;
%
%         function G=g(x)
%             G=tanh(x);
%
%         end
%         function G=gd(x)
%             G=1-(tanh(x)).^2;
%
%         end
% end
%     function y=initialspike(y,peak_interval)
%         Sp=zeros(size(y));label=ones(size(y,1),1);
%         for k=1:size(y,1)
%             %a=neo(y(k,:),1);
%             yp=y(k,:);yp(yp<0)=0;yn=-y(k,:);yn(yn<0)=0;
%             th1=spikethresh(yp);
%             th2=spikethresh(yn);
%             if th1>th2
%             Sp(k,:)=searchspike_in(yp,th1,peak_interval);
%             else
%             Sp(k,:)=searchspike_in(yn,th2,peak_interval);
%             end
%             a=diff(find(Sp(k,:)==1));
%             a((a>500|a<40))=[];
%             b=y(k,Sp(k,:)==1);
%             coef_s=std(a)/mean(a);
%             coef_w=std(b)/abs(mean(b));
%             if (coef_s<0.4)&&(coef_w<0.3)
%                 label(k)=0;
%             end
%
%         end
%                     y(label==1,:)=[];
%     end
% function [S,yt,xg,coef_s,coef_w,much]=ManualCICAtest(xw,y,peak_interval,fs,pic)
%         %              yp=y;yp(yp<0)=0;yn=-y;yn(yn<0)=0;
%         %             th1=spikethresh(yp);
%         %             th2=spikethresh(yn);
%         %             if th1>th2
%         %             S=searchspike_in(yp,th1,peak_interval);
%         %             else
%         %             S=searchspike_in(yn,th2,peak_interval);
%         %             end
%         if skewness(y')<0
%             y=-y;
%         end
%         th2=spikethresh(y);
%         S=searchspike_in(y,th2,peak_interval);
%         %%这个是debug用以后删掉
%         if isempty(S)
%             figure;plot(y);
%             title('bug');
%             pause
%         end
%         %%%%
%         [yt,xg]=fcica3(xw,S);
%         th2=spikethresh(yt);th1=0;c=1;
% %         if pic==1
%         figure;set(gcf,'Position',get(0,'ScreenSize'));
%         while abs(th1-th2)>1e-3 && c<20
%             th1=th2;
%             S=searchspike_in(yt,th1,peak_interval);
%             [yt,xg]=fcica3(xw,S);
%             %th2=spikethresh(yt);
%             [th2,lo]=spikethresh(yt);
%             plot(yt);hold on;plot(th2*ones(size(yt)));hold on;plot(lo*ones(size(yt)));pause(1);hold off;
%             c=c+1;
%         end
%         close
% %         else
% %             while abs(th1-th2)>1e-3 && c<20
% %             th1=th2;
% %             S=searchspike_in(yt,th1,peak_interval);
% %             [yt,xg]=fcica3(xw,S);th2=spikethresh(yt);
% %             c=c+1;
% %             end
% %         end
%
%
%         a=diff(find(S==1));
%         much=(length(a(a<1000))>=sum(a(a<1000))/fs*3);%发放个数的限制
%         a((a>500))=[];
%         b=yt(S==1);
%         %         if numel(b)>4
%         %         lillietest(b)
%         %         end
%         coef_s=std(a)/mean(a);
%         coef_w=std(b)/abs(mean(b));
%         %         xg
%         if c==20
%             disp('!!Warn:the automatic threshold does not converge');
%         end
%     end
end