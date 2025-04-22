% INVERSION SCRIPT %
clear;clc; % Clear command window

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% INPUTS %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== Data ======
sampleName = 'ICSV13';
model = 'inputs/forward_ICSV13.csv'; % Forward models.
measurements = 'inputs/input_ICSV13.csv'; % Measurements

% ====== Data type ======
dataFormat = 0; % What type of data do you have? 1 = all measurements. 0 = mean and std. of variables.
keep_cols = [1:13];% % These are the columns of the model, inlcuding T & P

% ====== Pressure units in forward model ======
unitsP = 1; % 1 = bar, Else = kbar.

% ====== Bootstrapping parameters ======
bootstrapType = 1;      % 1 = Parametric. Else = non-parametric.
it = 50;        % How many random iterations do you want to calculate?

% ====== Plotting ======
% Bootstrap progress
plotBoot = 0; % 1 = YES, else = NO.

% Inversion results
plotInv = 1; % 1 = YES, else = NO.
confidenceLevel = 0.68;  % Confidence level for 2D ellipse
boxplots = 1;   % Do you want boxplots or histograms? 1 = boxplot, 0 = histogram
plot_type = 0; % What type of plot do you want? 1 = contour plot, 0 = heatmap;
T_bins = 10; % Number of temperature bins in 2D histogram (Figure 2, 4)
P_bins = 10; % Number of pressure bins in 2D histogram (Figure 2, 4)

% Residuals
plotResiduals = 0; % 1 = YES, else = NO.

% Sensitivity
plotSens = 1; % 1 = YES, else = NO.



%%%%%%%%%%%%%%%%%%%%% CODE %%%%%%%%%%%%%%%%%%%%
%%%% BEST NOT TO ALTER UNLESS YOU ARE SURE %%%%
%%%%%%%%%%%%%%%%%%%%% CODE %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PART 1: Input the model and create the Pressure-Temperature grid %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in model
model = readtable(model,'VariableNamingRule','preserve'); model = sortrows(model,2);
model = model(:,keep_cols);

% Create Pressure-Temperature grid
data = table2array(model);
temperature = data(:,1);
pressure = data(:,2);
T = unique(temperature); T = T(~isnan(T)); ix = length(T);
P = unique(pressure); P = P(~isnan(P)); iy = length(P);
Tres = T(2) - T(1); 
Pres = P(2) - P(1);
variablesM = model.Properties.VariableNames;
[X,Y] = meshgrid(T,P);
model_data = data(:,3:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PART 2: Re-sample the observations %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check data compatibility
obs = readtable(measurements,'ReadRowNames',true,'VariableNamingRule','preserve');
rows = obs.Properties.RowNames;
if ~any(strcmp(rows, 'MEAN')) && dataFormat == 0
    error('This data does not contain MEAN or STD information. Wrong data type.')
elseif any(strcmp(rows, 'MEAN')) && dataFormat == 1
    error('This data contains MEAN or STD information. Wrong data type.')
end

% Read in data
if dataFormat == 0
    obs = readtable(measurements,'ReadRowNames',true,'VariableNamingRule','preserve');
    rows = obs.Properties.RowNames;
else
    obs = readtable(measurements,'VariableNamingRule','preserve');
end
obs = obs(:,keep_cols(3:end)-2);
variables = obs.Properties.VariableNames;
obs = table2array(obs);

% Check whether columns of forward model match columns of measurements
if ~isequal(variablesM(3:end),variables)
    error(['The column names for the model and the measurements are different. ' ...
        'Check the model is correct.'])
else
    clear variablesM
    variables = cellfun(@(x) strrep(x, '_', ''), variables, 'UniformOutput', false);
end

% Find mean and standard deviation of measurements
if dataFormat == 0
    sigma = obs(2,:); mu = obs(1,:);
else
    sigma = std(obs,1); mu = mean(obs,1);
end

% Bootstrap re-sampling
if bootstrapType == 1 % Parametric
    samples = Functions_NO_EDIT.gaussian_boot(it,mu,sigma);
else % Non-parametric
    samples = Functions_NO_EDIT.nonpara_boot(it,obs);
end
samples = [samples;mean(samples,1)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PART3: Perform the grid-search inversion and sensitivity analysis %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid-search inversion
[t_best,p_best,model_misfit] = Functions_NO_EDIT.gridSearch(it,samples,model_data,temperature,pressure);

% Find mean and std solution
tMean = mean(t_best); pMean = mean(p_best);
stdT = std(t_best); stdP = std(p_best);

% Find median and IQR solution
tMed = median(t_best); pMed = median(p_best);
iqrT1 = prctile(t_best,25); iqrT2 = prctile(t_best,75);
iqrP1 = prctile(p_best,25); iqrP2 = prctile(p_best,75);

% Find fit metrics
distances = sqrt((temperature - tMed).^2 + (pressure - pMed).^2);
[~, closestRow] = min(distances);
model_prediction = model_data(closestRow,:);
d = mu - model_prediction;
chi_variables = abs(d) ./ (2 * sigma); chiMedVar = chi_variables;
chiMed = sum(abs(d) ./ (2 * sigma))/length(d);

% Compute sensitivity
[T_max,T_min,P_max,P_min] = Functions_NO_EDIT.sensitivity(model_data,samples,tMean,pMean,temperature,pressure,sigma,mu,bootstrapType,dataFormat,it);
tmp01 = [abs(T_max);abs(T_min)]; tmax = max(tmp01,[],1);
tmp01 = [abs(P_max);abs(P_min)]; pmax = max(tmp01,[],1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PART4: Print diagnostics and plot results %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correct units
if unitsP == 1
    p_best = p_best./1000;
    Y = Y./1000;
    P = P./1000;
    P_min = P_min/1000;
    P_max = P_max/1000;
    pMean = pMean/1000;
    stdP = stdP/1000;
    pMed = pMed/1000;
    iqrP1 = iqrP1/1000;
    iqrP2 = iqrP2/1000;
    Pres = Pres/1000;
    pmax = pmax./1000;
end

% Variables fit
results = [keep_cols(3:end); mu + 2*sigma; mu - 2*sigma; sigma./mu *100; model_prediction; chiMedVar; tmax; pmax];
dfit = array2table(results,'VariableNames',variables,'RowNames',{'Column','μ + 2σ', 'μ - 2σ','σ % of μ','Prediction','X_i','ΔT (°C)','ΔP (kbar)'});

% Save results
filename = "output_variables/TPsolutions_" + sampleName + ".csv";
writematrix(["Temperature", "Pressure"; [t_best(:), p_best(:)]], filename);


% Plot results
if plotBoot == 1
    Functions_NO_EDIT.plotBoot(t_best,p_best,sampleName)
end
if plotInv == 1
    Functions_NO_EDIT.plotInv(t_best,p_best,X,Y,P,T,model_misfit,confidenceLevel,ix,iy,plot_type,boxplots,T_bins,P_bins,sampleName)
end
if plotResiduals == 1
    Functions_NO_EDIT.plotResiduals(mu,sigma,model_prediction,obs,dataFormat,variables,sampleName)
end
if plotSens == 1
    Functions_NO_EDIT.plotSens(T_min,T_max,P_min,P_max,variables,tMed,pMed,sampleName)
end

% Write inversion diagnostics
no_variables = width(obs);
fit_var = length(find(chiMedVar < 1));

score(1) = chiMed;
score(2) = stdT/sqrt(no_variables);
score(3) = stdP/sqrt(no_variables);
score(4) = no_variables;
score(5) = fit_var;
score(6) = Tres;
score(7) = Pres;
score(8) = it;
clc;
Functions_NO_EDIT.writeResults(tMean,stdT,pMean, stdP, tMed,iqrT1,iqrT2,pMed,iqrP1,iqrP2,score,dfit,sampleName)

%%%%%%%%%%%%%%%%%%%%%
%%%% END OF CODE %%%%
%%%%%%%%%%%%%%%%%%%%%
