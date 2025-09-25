function [mis,dT,dP,distT,distP,Tbest,Pbest] = LOOfit(model_data,samples,temperature,pressure,tMed,pMed,mu,sigma)


nCols = size(model_data,2);
nRows = size(model_data,1);
errors = zeros(1, nCols);  % prediction error for left-out column
distT = zeros(1, nCols);
distP = zeros(1, nCols);
Tbest = zeros(1, nCols);
Pbest = zeros(1, nCols);
for i = 1:nCols
    % Indices excluding i
    idx = setdiff(1:nCols, i);

    % Model subset and data subset
    modelData = model_data(:, idx);
    data = samples(end, idx);

    % Find row with smallest misfit (L1 norm) on training subset
    normalizedDiff = abs(modelData - data) ./ abs(data);
    misfit = sum(normalizedDiff, 2);

    % Find row with smallest normalized misfit
    [v, rowIdx] = min(misfit);

    % Best (T, P) from training subset
    bestT = temperature(rowIdx);
    bestP = pressure(rowIdx);

    % --- Predict left-out column ---
    predictedVal = model_data(rowIdx, i);
    trueVal      = samples(end, i);

    % Prediction error for left-out column
    errors(i) = abs(predictedVal - trueVal)./trueVal*100;
    dpred = model_data(rowIdx,idx);
    dmeasured = mu(idx);
    sigma_id = sigma(idx);
    d = dpred - dmeasured;
    chiV = sum(abs(d) ./ (2 * sigma_id))/length(d);
    mis(i) = chiV;

    % Store T, P as well
    dT(i) = bestT;
    dP(i) = bestP;
    distT(i) = bestT - tMed;
    distP(i) = bestP - pMed;
    Tbest(i) = bestT;
    Pbest(i) = bestP;
end