clc
clear
close all
load('1y.mat')
[cor1 lag1]=xcorr(y,'unbiased');
figure(1)
plot(lag1,cor1)
sum(y.*y)