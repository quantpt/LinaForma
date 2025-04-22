% This script will produce tornado sensitivity plots for each parameter
clear;clc;

%%%%%%%%% INPUTS %%%%%%%%%
% ====== Data ======
model = 'inputs/forward_model.csv';
measurements =  'inputs/InputA.csv';

% ====== Data type ======
raw = 0; % What type of data do you have? 1 = all measurements. 0 = mean and std. of variables.

% ====== Bootstrapping parameters ======
bootstrap_type = 1;      % What type of bootstrapping do you want to use? Parametric = 1, non-parametric = 0
it = 500;        % How many random iterations do you want to calculate? The more then better.

% ====== Select P-T point for sensitivity analysis ======
T_best = 570;      % Select the best-fit T point to which results will be compared.
P_best = 9300;      % Select the best-fit P point to which results will be compared.


%%%%%%%%%%%%%%%%%%%%% CODE %%%%%%%%%%%%%%%%%%%%
%%%% BEST NOT TO ALTER UNLESS YOU ARE SURE %%%%

%%%% PART 1: Input the data and create the Pressure-temperature grid %%%%
% Read in data
model_tmp = readtable(model); model_tmp = sortrows(model_tmp,2);

% Create Pressure-Temperature grid
data = table2array(model_tmp);
temperature = data(:,1);
pressure = data(:,2);
T = unique(temperature); T = T(~isnan(T)); ix = length(T);
P = unique(pressure); P = P(~isnan(P)); iy = length(P);
[X,Y] = meshgrid(T,P);
model_data = data(:,3:end);

%%%% PART 2: Boostrapping the observations, keeping all variables but 1 constant %%%%
if raw == 1
    obs = readtable(measurements);
    var_tmp = obs.Properties.VariableNames;
    if any(strcmp(var_tmp, 'Statistic'))
        error('This data only contains MEAN or STD information.')
    end
    obs = table2array(obs);
elseif raw == 0
    obs = readtable(measurements,'ReadRowNames',true);
    var_tmp = obs.Properties.RowNames;
    if ~any(strcmp(var_tmp, 'MEAN'))
        error('This data does not contain MEAN or STD information.')
    end
    syn_mean = table2array(obs(1,:)); syn_sigma =table2array(obs(2,:));
    obs = table2array(obs);
end


% Loop through each variable
for n_variable = 1:size(model_data,2)

% Perform bootstrap re-sampling
if bootstrap_type == 1 % Parametric bootstrapping
    if raw == 0
        sigma = syn_sigma; mu = syn_mean;
    else
        sigma = std(obs,1); mu = mean(obs,1);
    end
    samples1 = +Functions_NO_EDIT.gaussian_boot(it,mu,sigma);
    samples = repmat(mu,it,1);
    samples(:,n_variable) = samples1(:,n_variable);

elseif bootstrap_type == 0 % Non-parametric
    if raw == 0
        error('You cannot perform non-parametric bootstrapping on this dataset.')
    end
    samples1 = bootstrp(it,@mean,obs);
    mean_samples = mean(samples1,1);
    samples = repmat(mean_samples,it,1);
    samples(:,n_variable) = samples1(:,n_variable);
end

%%%% PART 3: Perform the grid-search inversion for each bootstrap resample %%%%
loop1 = 0;
for j = 1:size(samples,1)
loop1 = loop1 + 1; sample_values = samples(j,:);
for ii = 1:length(temperature)
    model_values = model_data(ii,:); 
    model_misfit(ii,:) = +Functions_NO_EDIT.misfit(model_values,sample_values);  
end
t_best(loop1,:) = temperature(model_misfit == min(model_misfit)); 
p_best(loop1,:) = pressure(model_misfit == min(model_misfit));
end

% Find variability of each variable
T_max(n_variable) = max(t_best)-T_best;
T_min(n_variable) = T_best - min(t_best);
P_max(n_variable) = max(p_best) - P_best;
P_min(n_variable) = P_best - min(p_best);

end

%%%% PART 4: Plots %%%%
variables = readtable(model); variables = sortrows(variables,2); variables = variables.Properties.VariableNames;
variables = variables(:,3:end);

% Plot T variability
fig1 = figure(1);
set(fig1,'Units','centimeters')
set(fig1,'Position',[0 0 0.9*21 0.3*29.7])
subplot(1,2,1)
impact = [-T_min;T_max]'; % Impact on the outcome (positive or negative)

% Sort variables based on their impact
positive_rows = all(impact > 0, 2);
impact(positive_rows, 1) = 0;
negative_rows = all(impact < 0, 2);
impact(negative_rows, 2) = 0;
total_impact = sqrt((impact(:,1) - impact(:,2)).^2);

[~, sortOrder] = sort(total_impact, 'descend');
sortedVariables = variables(sortOrder);
sorted_impact = impact(flipud(sortOrder),:);

% Create tornado plot
barh(sorted_impact,'stacked'); % Blue bars for positive impact
hold on;
set(gca, 'YTick', 1:length(variables), 'YTickLabel', sortedVariables,'Fontsize',12); % Set y-axis labels
xlabel('Temperature (°C)'); % Label x-axis
t = append('Variation relative to ',string(T_best),' °C');
title(t); % Title
xlim([min(min(sorted_impact))-10 max(max(sorted_impact))+10])


subplot(1,2,2)
impact = [-P_min;P_max]'; % Impact on the outcome (positive or negative)

% Sort variables based on their impact
positive_rows = all(impact > 0, 2);
impact(positive_rows, 1) = 0;
negative_rows = all(impact < 0, 2);
impact(negative_rows, 2) = 0;
total_impact = sqrt((impact(:,1) - impact(:,2)).^2);

[sortedImpact, sortOrder] = sort(total_impact, 'descend');
sortedVariables = variables(sortOrder);
sorted_impact = impact(flipud(sortOrder),:);

% Create tornado plot
barh(sorted_impact,'stacked');
set(gca, 'YTick', 1:length(variables), 'YTickLabel', sortedVariables,'Fontsize',12); % Set y-axis labels
xlabel('Pressure (bar)'); % Label x-axis
t = append('Variation relative to ',string(P_best),' bar');
title(t)
lg = legend('-','+','Location','southeast');
xlim([min(min(sorted_impact))-100 max(max(sorted_impact))+100])

% Save figure
saveas(fig1,"FIGURES/L4_fig1.pdf");
disp('FINISHED')










