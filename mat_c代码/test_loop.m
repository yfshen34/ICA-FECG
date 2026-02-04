function [sumgd2]=test_loop(y,k)
    sumgd2=0;
     for k=1:2999
        tkt = circshift(y',-k);
        tkt(end-k+1:end)=0;
        tk = circshift(y',k);
        tk(1:k)=0;
        T=tk+tkt;
        sumgd2=sumgd2+mean(T);
     end
 return 