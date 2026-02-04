function lambda = genCoolingFF(t,gamma,lambda_0)
    % lambda = lambda_0 / sample^gamma
    lambda = lambda_0 ./ (t .^ gamma);
end