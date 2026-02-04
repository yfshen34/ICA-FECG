%% 评估abdominal-and-direct-fetal-ecg-database-1.0.0
function [acc,ppv,sen,f1]=evaluation(S,Findex,selection,t)
% t:评估序列的时间
% %S:1*sample
% %Findex:ann
%     t=1:size(record_filter,2);
%     S_index=S;
% load('all_r01.mat');
Fs=1000;
if nargin<4
    t=60000;
end
if selection==1
    Sindex=S;
else
    Sindex=find(S==1);
end

F=zeros(1,t);
for i=1:size(Findex,1)
    j=single(Findex(i));
    F(j)=1;
%     Findex=[Findex j];
end
win=Fs/200;
Findex(Findex<win+1)=[];
Sindex(Sindex<win+1)=[];
Findex(Findex>t-win+1)=[];
Sindex(Sindex>t-win+1)=[];
length(Findex);
% 匹配
TP=0;%correctly detected FECG QRS complexes
for i=1:length(Findex)
    temp=abs(Sindex-Findex(i));
    [num,loc]=min(temp);
    if isempty(num)
        aa=1;
    end
    if ~isempty(num) && num(1)<=50% 阈值<50ms
        Sindex(loc(1))=[];% 消去
        Findex(i)=-1;
        TP=TP+1;
    end
end
FP=length(Sindex);%falsely detected nonexistent FECG QRS complexes
FN=length(find(Findex~=-1));%falsely missed existent FECG QRS complexes


% acc
acc=100*TP/(TP+FP+FN);
% ppv
ppv=100*TP/(TP+FP);
% sen
sen=100*TP/(TP+FN);
% f1
f1=(2*ppv*sen)/(ppv+sen);

% plot(S)
% hold on
% plot(Findex,1,'r+')

%FHR 每分钟的心率
% for i=1:size(Findex,2)/(60*Fs)
%     seg1=S(2,(i-1)*60*Fs+1:i*60*Fs);
%     fhr(i)=length(find(seg1==1));
%     seg2=F((i-1)*60*Fs+1:i*60*Fs);
%     ref(i)=length(find(seg2==1));
% end

% string=strcat('eva_',s{1},'.mat');
% save(string, 'sen','acc','ppv','f1','fhr','ref');

