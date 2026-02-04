function state = dynamicWhitening(blockdata, dataRange, state, adaptiveFF)

    nPts = size(blockdata,2);

    % define adaptive forgetting rate: lambda
    switch adaptiveFF.profile
        case 'cooling'
            lambda = genCoolingFF(state.counter+dataRange, adaptiveFF.gamma, adaptiveFF.lambda_0);
            if lambda(1) < adaptiveFF.lambda_const
                lambda = repmat(adaptiveFF.lambda_const,1,nPts); 
            end
        case 'constant'
            lambda = repmat(adaptiveFF.lambda_const,1,nPts);
        case 'adaptive'
            lambda = repmat(state.lambda_k(end),1,nPts); % using previous adaptive lambda_k from adaptiveOrica
    end
        
    % update sphere matrix using online RLS whitening block update rule
    v = state.icasphere * blockdata; % pre-whitened data 
    lambda_avg = 1 - lambda(ceil(end/2));    % median lambda
    QWhite = lambda_avg/(1-lambda_avg) + trace(v' * v) / nPts;
    state.icasphere = 1/lambda_avg * (state.icasphere - v * v' / nPts / QWhite * state.icasphere);

end
