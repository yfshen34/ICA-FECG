function lambda = genAdaptiveFF(dataRange,lambda,decayRateAlpha,upperBoundBeta,transBandWidthGamma,transBandCenter,ratioOfNormRn)
% lambda = lambda - DecayRate*lambda + UpperBound*Gain*lambda^2
% Gain(z) ~ tanh((z/z_{min} - TransBandCenter) / TransBandWidth)
    gainForErrors = upperBoundBeta*0.5*(1+tanh((ratioOfNormRn-transBandCenter)/transBandWidthGamma));
    f = @(n) (1+gainForErrors).^n * lambda(end) - decayRateAlpha*((1+gainForErrors).^(2*n-1)-(1+gainForErrors).^(n-1))/gainForErrors*lambda(end)^2;
    lambda = f(1:length(dataRange));    
end
