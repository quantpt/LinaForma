% This script plots the residuals for a particular point of interest
clear;clc;

%%%%%%%%% INPUTS %%%%%%%%%
% ====== Data ======
model = 'inputs/forward_model.csv'; % Forward models.
measurements = 'inputs/InputA.csv'; % Measurements

% ====== Data type ======
raw = 0; % What type of data do you have? 1 = all measurements. 0 = mean and std. of variables.

% ====== Select P-T point for forward model data ======
T_best = 590; % Select the T point to which results will be compared.
P_best = 9300;      % Select the P point to which results will be compared. Units = bars.



%%%%%%%%%%%%%%%%%%%%% CODE %%%%%%%%%%%%%%%%%%%%
%%%% BEST NOT TO ALTER UNLESS YOU ARE SURE %%%%
model_tmp = readmatrix(model); model_tmp = sortrows(model_tmp,2);
if raw == 1
    measurements_tmp = readtable(measurements);
    var_tmp = measurements_tmp.Properties.VariableNames;
    if any(strcmp(var_tmp, 'Statistic'))
        error('This data only contains MEAN or STD information.')
    end
    measurements = readmatrix(measurements);
    tmp = size(measurements,2);
elseif raw == 0
    syn = readtable(measurements,'ReadRowNames',true);
    var_tmp = syn.Properties.RowNames;
    if ~any(strcmp(var_tmp, 'MEAN'))
        error('This data does not contain MEAN or STD information.')
    end
    syn_mean = table2array(syn(1,:)); syn_sigma =table2array(syn(2,:));
    tmp = size(syn_mean,2);
end
variables = readtable(model); variables = sortrows(variables,2); variables = variables.Properties.VariableNames;
variables = variables(:,3:end);


% Get data
m_data = model_tmp(:,3:end); TP = model_tmp(:,1:2); % Get PT points of model
row = +Functions_NO_EDIT.find_points([T_best,P_best],TP,0);
data = m_data(row,:);


for i = 1:tmp

        if raw == 0
            mu = syn_mean(i); sig = syn_sigma(i);
            x = linspace(mu-3*sig, mu + 3*sig, 1000); 
            distribution = normpdf(x, mu, sig);
            sigma(i) = sig; mn(i) = mu;
            observed_dist(i,:) = distribution;
            x_axis(i,:) = x;

        else
            tmp = measurements(:,i);
            sigma(i) = std(tmp);
            mn(i) = mean(tmp);
            [f,x] = ksdensity(tmp);
            observed_dist(i,:) = f;
            x_axis(i,:) = x;
        end

end

fig1 = figure(1);
row = ceil(length(variables)/3);
t = tiledlayout('flow');
set(fig1,'Units','centimeters')
set(fig1,'Position',[0 0 0.9*21 row*1/5*29])
for i = 1:length(variables)
    nexttile
    plot(x_axis(i,:),observed_dist(i,:)); hold on
    plot([data(1,i),data(1,i)],[0 max(observed_dist(i,:))],'r--','LineWidth',2)
    if raw ~= 0
        plot([min(measurements(:,i)) min(measurements(:,i))],[0 max(observed_dist(i,:))])
        plot([max(measurements(:,i)) max(measurements(:,i))],[0 max(observed_dist(i,:))])
    else
        plot([(syn_mean(:,i)-2*syn_sigma(:,i)) (syn_mean(:,i)-2*syn_sigma(:,i))],[0 max(observed_dist(i,:))])
        plot([(syn_mean(:,i)+2*syn_sigma(:,i)) (syn_mean(:,i)+2*syn_sigma(:,i))],[0 max(observed_dist(i,:))])
    end
    mnT = round(mn(i),3,'significant');
    sigT = round(sigma(i),2,'significant');
    t = append('Mean = ',string(mnT),' ± ',string(sigT),' (1σ)');
    title(t);
    ylabel('P.D.E.')
    xlabel(string(variables(i)))
end
% if raw == 0
%     lg = legend({'Observations','Model result','μ - 2σ','μ + 2σ'});
% elseif raw ~= 0
%     lg = legend({'Observations','Model result','Minimum observation','Maximum observation'});
% end
% lg.Layout.Tile = 'South';
print('-vector',fig1,"FIGURES/L3_fig1.pdf",'-dpdf','-bestfit')
disp('FINISHED')
