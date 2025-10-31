% INVERSION SCRIPT %
clear;clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% INPUTS %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== Data ======
sampleName = 'ICSV13';
model = 'inputs/tableS2_forward_ICSV13.csv'; % Forward models.
measurements = 'inputs/tableS1_measurements_ICSV13.csv'; % Measurements

% ====== Data type ======
dataFormat = 0; % What type of data do you have? 1 = all measurements. 0 = mean and std. of variables.
keep_cols = [1,2,3,4,5,6,7,8,9,10,11]; % These are the columns of the model, inlcuding T & P

% ====== Pressure units in forward model ======
unitsP = 1; % 1 = bar, Else = kbar.

% ====== Bootstrapping parameters ======
bootstrapType = 1;      % 1 = Parametric. Else = non-parametric.
it = 50;        % How many random iterations do you want to calculate?

% ====== Plotting ======
% Bootstrap progress
plotBoot = 0; % 1 = YES, else = NO.

% Inversion results
plotInv = 0; % 1 = YES, else = NO.
confidenceLevel = 0.68;  % Confidence level for 2D ellipse
boxplots = 0;   % Do you want boxplots or histograms? 1 = boxplot, 0 = histogram
plot_type = 1; % What type of plot do you want? 1 = contour plot, 0 = heatmap;
T_bins = 10; % Number of temperature bins in 2D histogram (Figure 2, 4)
P_bins = 10; % Number of pressure bins in 2D histogram (Figure 2, 4)

% Residuals
plotResiduals = 0; % 1 = YES, else = NO.

% Sensitivity
plotSens = 1; % 1 = YES, else = NO.

% Ftotal(i)
plotLOOA = 0; % 1  = YES, else = NO.



%%% CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Functions_NO_EDIT.run_script(sampleName, model, measurements, dataFormat, keep_cols, unitsP, ...
                   bootstrapType, it, plotBoot, plotInv, confidenceLevel, ...
                   boxplots, plot_type, T_bins, P_bins, plotResiduals, plotSens, plotLOOA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

