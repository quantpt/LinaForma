function [T_max,T_min,P_max,P_min] = sensitivity(model_data,obs,T_best,P_best,temperature,pressure,syn_sigma,syn_mean,bootstrap_type,raw,it)

% Loop through each variable
for n_variable = 1:size(model_data,2)

% Perform bootstrap re-sampling
if bootstrap_type == 1 % Parametric bootstrapping
    if raw == 0
        sigma = syn_sigma; mu = syn_mean;
    else
        sigma = std(obs,1); mu = mean(obs,1);
    end
    samples1 = Functions_NO_EDIT.gaussian_boot(it,mu,sigma);
    samples = repmat(mu,it,1);
    samples(:,n_variable) = samples1(:,n_variable);

else % Non-parametric
    if raw == 0
        error('You cannot perform non-parametric bootstrapping on this dataset.')
    end
    samples1 = bootstrp(it,@mean,obs);
    mean_samples = mean(samples1,1);
    samples = repmat(mean_samples,it,1);
    samples(:,n_variable) = samples1(:,n_variable);
end

%%%% PART 3: Perform the grid-search inversion for each bootstrap resample %%%%
[t_best,p_best,~] = Functions_NO_EDIT.gridSearch(it,samples,model_data,temperature,pressure);

% FInd mean and standard dev.
tmean = mean(t_best); pmean = mean(p_best);
stdT = std(t_best); stdP = std(p_best);
tMAX = tmean + 2*stdT;
tMIN = tmean - 2*stdT;
pMAX = pmean + 2*stdP;
pMIN = pmean - 2*stdP;


% Find variability of each variable
T_max(n_variable) = tMAX-T_best;
T_min(n_variable) = T_best - tMIN;
P_max(n_variable) = pMAX - P_best;
P_min(n_variable) = P_best - pMIN;

end