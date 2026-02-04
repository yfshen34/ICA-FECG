function [sign]=preprocess(s)
    x=s(:,2:end);
    y=s(:,1:end-1);
    [A,B,r,U,V,stats]=canoncorr(x',y');
    r=r';
     multiPlot('CCA',U')
%     decide(U,A,r,x);
        th=0.97;
        AA=U;
        for com=1:4
            if r(com)<th
                AA(:,com)=0;
            end
        end
%         ÖØ¹¹ÐÅºÅX=A'*x
        sign=AA*inv(A);
        sign=sign';
        multiPlot('CCA',sign)
end

