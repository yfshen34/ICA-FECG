function S=searchspike_in(x,high,peak_interval,fs,graph,ann,outlier)
%————找到x中峰值的位置————
if nargin<4
    fs=1000;
end

if nargin<5
    graph=0;
end

if nargin<7
    outlier=15;
end
win=0.005*fs;
S=zeros(size(x));
win=fix(fs/200);
x_jd=x(2*win+1:end-2*win);

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
    plot(high*ones(1,length(X)),'r-');
%     legend('选取constrainted FastICA的spike去趋势后：');
    hold off
end

if high==0
    S=[];
    return
elseif high>0
%     yp=x_jd;
    yp=abs(x_jd);
else
%     yp=-x_jd;
    yp=abs(x_jd);
end


% 先找异常值，若有，看是否附近有超过阈值的Spike
[~,L2] = findpeaks(yp,'minpeakheight',outlier);%x的rms=1
[~,L] = findpeaks(yp,'minpeakheight',high,'MINPEAKDISTANCE',peak_interval);

if ~isempty(L2)
    % 有超过阈值的Spike，用临近点代替
    [value,L3] = findpeaks(yp,'minpeakheight',high);
    for k=1:length(L2)
        temp=abs(L2(k)-L3);
        [num,loc]=sort(temp,'ascend');%距离从小到大
        i=1;
        while num(i)<peak_interval/2% 阈值<peak_interval且不为outlier
        if value(loc(i))<outlier
            L(L==L2(k))=L3(loc(1));
            break
        else
            i=i+1;
        end
        end
    end
end

S(1,L)=1;

if graph
    hold on
    plot(L+2*win,ones(1,length(L+2*win)),'k+');
    
    
    
    if nargin==6
        ann1=ann-2*win;
        plot(ann1,-1,'r+');
        figure
        plot(S)
        hold on
        plot(ann,1,'r+');
        legend('根据阈值选取的spike');
    end
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
%             S(1,L3)=2;
% %             S(1,L)=1;
%         else S=[];
%         end
end