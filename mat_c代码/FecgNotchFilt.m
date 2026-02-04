function Xf=FecgNotchFilt(X,fs,cName,graph,dbFlag)
% -----------------------------------------------------------------------------------
%   ECG power line canceling by notch filtering
%
% Xf=FecgNotchFilt(X,fs,cName,graph,dbFlag)
%
% X          : input signal matrix (one signal per column)
% fs         : sampling frequency
% cName      : record name
% graph      : flag enabling figure drawing
% Xf         : 
%
% Author: Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
% For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% -----------------------------------------------------------------------------------

if(nargin<3), cName=''; end
if(nargin<4), graph=0; end
if(nargin<5), dbFlag=0; end
graphD= graph &dbFlag ;
graphSpt= graph &dbFlag ;
%-------------------------------------------------------------
% recording duration
[ndt, ns]=size(X);
vtime= [1:ndt]/fs;
Xf=zeros(size(X));
for is=1:ns
    if(dbFlag && graphSpt), figure; pwelch(X(:,is),[],[],[],fs);
        title('Detrended ECG Welch spectrum');
    end
    [Px,Fv]=pwelch(X(:,is),[],[],[],fs);
    iF50w=find(49<Fv & Fv<51);  F50w=Fv(iF50w); % 检测50Hz干扰
    [maxP50,imaxP50]=max(Px(iF50w));            %找到50Hz频段的最大功率点
    f50mP=F50w(imaxP50);                        %实际干扰频率（可能非精确50Hz)
    iFes50w=find((46<Fv & Fv<49) | (51<Fv & Fv<54));
    pem50w=mean(Px(iFes50w));                   %领域平均功率
    pestd50w=std(Px(iFes50w));                  %领域功率标准差
    d50P=maxP50-pem50w;                         %50Hz功率与领域的差异
    e50=d50P>5*pestd50w;                        %差异>5倍标准差
    
    iF60w=find(59<Fv & Fv<61);  F60w=Fv(iF60w); %同理检测60Hz
    [maxP60,imaxP60]=max(Px(iF60w));
    f60mP=F60w(imaxP60);
    iFes60w=find((56<Fv & Fv<59) | (61<Fv & Fv<64));
    pem60w=mean(Px(iFes60w)); pestd60w=std(Px(iFes60w));
    d60P=maxP60-pem60w;
    e60=d60P>5*pestd60w;
    
   % if(e50 || e60)
    if any(e50) || any(e60)
        if(d50P>d60P), fnotch=f50mP;
        else fnotch=f60mP;  end
        
        % Notch filter
        fnotchn=fnotch/fs; bwn=1/fs; grafSptn=0;
        
        Xf(:,is)= filtNotchFB(X(:,is),fnotchn,bwn, grafSptn);
        Xf(:,is)= filtNotchFB(Xf(:,is),2*fnotchn,bwn, grafSptn);
        Xf(:,is)= filtNotchFB(Xf(:,is),3*fnotchn,bwn, grafSptn);
        Xf(:,is)= filtNotchFB(Xf(:,is),4*fnotchn,bwn, grafSptn);
    else
        Xf(:,is)=X(:,is);
    end
end
if(graph)
    figure, set(gcf,'Color','white');
    for is=1:ns,
        subplot(ns,1,is), hold on, plot(vtime,X(:,is),'r'), plot(vtime,Xf(:,is),'b');
        wgmimaV=mimaxscG(X(:,is),0,0,.1);
        ylim(wgmimaV);
        if(is~=ns), set(gca,'XTickLabel',''); end
        if(is==1), title([cName,': original(r) & notch filtered(b) ECG']); end
    end
    shg
end

return
end %== function ================================================================
%


