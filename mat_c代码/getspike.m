function [F,x,thres,coef_s,coef_w,ind,Spike,Sc,L,C]=getspike(x,fs,peakinterval,graph,ann,th)%解决符号不一致的一些问题
%         interval=60/300*fs;%fix是靠0取整
%         interval=fix(fs/100);%fix是靠0取整
% 去头去尾
win=fix(1000/100);
% x_jd=x(2*win+1:end-2*win);
% x_jd=x;

% 归一化
% x_jd=x_jd/prctile(abs(x_jd),99);
% x_sign=x_jd;
% x_jd=power(x_jd,2);
% x_sign(find(x_sign>0))=1;
% x_sign(find(x_sign<0))=-1;
% x_jd=x_jd.*x_sign;

% 去趋势
% in=x_jd';
% fin=[in(1:fs);in];
% np = length(fin);
% xd = decimate(fin,round(fs/50),'fir'); % decimation
% lbx = medfilt1(xd,10); % median filtering - get isoline
% lb = interp(lbx,round(fs/50)); % interpolation
% out = fin-lb(1:np);
% lbx1=medfilt1(fin,500);
% out1=fin-lbx1(1:np);
% out(1:fs)=[];  % result
% X=out';
X=x';


yp=X(2*win+1:end-2*win);
yn=-X(2*win+1:end-2*win);


if nargin==6
    if th>0
    [~,L1] = findpeaks(yp,'minpeakheight',th,'MINPEAKDISTANCE',peakinterval);
    [~,L2] = findpeaks(yp,'minpeakheight',18,'MINPEAKDISTANCE',peakinterval);%x的rms=1
    L3=setdiff(L1,L2);%返回在L1中，不在L2中的元素
    th1=sum(x_jd(L3));
    L=L3;
    
    else
    [~,L1i] = findpeaks(yn,'minpeakheight',-th,'MINPEAKDISTANCE',peakinterval);
    [~,L2i] = findpeaks(yn,'minpeakheight',18,'MINPEAKDISTANCE',peakinterval);%x的rms=1
    L3i=setdiff(L1i,L2i);%返回在L1中，不在L2中的元素
    th2=sum(x_jd(L3i))
    x=-x;
    L=L3i;
    
    end
    
else
    yp_sort=sort(yp);
    th1=spikethresh(yp_sort(1:fix(0.995*length(yp_sort))));
    yn_sort=sort(yn);
    th2=spikethresh(yn_sort(1:fix(0.995*length(yn_sort))));
    %         th1=spikethresh(yp_fil)
    %         th2=spikethresh(yn_fil)
    
    if th1>th2
        %通过分位点确定符号
        %     prctile(yp,90)
        %     prctile(yn,90)
        %     if prctile(yp(L3),90)>prctile(yn(L3i),90)
        
        % S=searchspike_in(yp,th1,peak_interval);
        [~,L1] = findpeaks(yp,'minpeakheight',th1,'MINPEAKDISTANCE',peakinterval);
        %             [~,L1] = findpeaks(yp,'minpeakheight',th1,'MINPEAKDISTANCE',peakinterval);
        [~,L2] = findpeaks(yp,'minpeakheight',18,'MINPEAKDISTANCE',peakinterval);%x的rms=1
        L3=setdiff(L1,L2);%返回在L1中，不在L2中的元素
        %             th1=mean(x_jd(L3));
        L=L3;
        thres=th1;
    else
        
        % S=searchspike_in(yn,th2,peak_interval);
        [~,L1i] = findpeaks(yn,'minpeakheight',th2,'MINPEAKDISTANCE',peakinterval);
        %             [~,L1i] = findpeaks(yn,'minpeakheight',th2,'MINPEAKDISTANCE',peakinterval);
        [~,L2i] = findpeaks(yn,'minpeakheight',18,'MINPEAKDISTANCE',peakinterval);%x的rms=1
        L3i=setdiff(L1i,L2i);%返回在L1中，不在L2中的元素
        %             th2=mean(x_jd(L3i));
        L=L3i;
        x=-x;
        X=-X;
        thres=th2;
    end
    
end
% 
% if graph
%     hold on;
%     plot(thres*ones(1,length(x)),'r-');
%     %  yline(th,'r--');hold on;
%     plot(L+2*win,1,'g+')
% end

% 恢复
L=L+2*win;

% spike小于下限
if numel(L)<length(x(1,:))/fs/1.5  %返回L中元素的个数
    F=[];F2=[];Spike=[];Sc=[];L=[];C=[];coef_s=[];coef_w=[];ind=0;
    return;
end

%140ms,100ms
left=win*3;
right=win*3;
Spike=zeros(length(L),left+right+1);
for ii=1:length(L)
%     [~,pos]=max(X(L(ii)-win:L(ii)+win));
%     pos=pos+L(ii)-win-1;
    pos=L(ii);
    if L(ii)-left<=0
        tt=X(1:pos+right);
        Spike(ii,left+right+2-length(tt):end)=tt;
    elseif L(ii)+right>length(x)
        tt=X(pos-left:end);
        Spike(ii,1:length(tt))=tt;
    else
        
    Spike(ii,:)=x(pos-left:pos+right);
    end
    
    
end
Sc(:,1)=abs(max(Spike,[],2)-min(Spike,[],2));%高
Sc(:,2)=sqrt(mean(Spike.^2,2));%每一行的平方根
% Sc(:,2)=sqrt(mean(diff(Spike,1,2).^2,2));%diff差分运算，1，2指每行后面的减去前面的
%             figure;
%             scatter(Sc(:,1),Sc(:,2));
% [C]=valley_seeking(Sc,17,0.3,3);
% 
% 
% F=[];
% F2=cell(1,max(C));
index=1;
% for ii=1:max(C)
%     F2{ii}=L(C==ii);
%     if length(L(C==ii))>length(x(1,:))/fs/1.5
        F{index}=L;
        a=diff(F{index});
        % much=(length(a(a<1000))>=sum(a(a<1000))/fs*3);%发放个数的限制
        a((a>500))=[];
        b=x(index,F{index});
        coef_s(index)=std(a)/mean(a);
        coef_w(index)=std(b)/abs(mean(b));
        index=index+1;
%     end
% end
ind=index-1;

end

function y=normalize(x)
s=max(abs(x));
y=x/s;
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
