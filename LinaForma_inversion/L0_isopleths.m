%% This script will show overlap between variables in P-T space
clear;clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% INPUTS %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== Data ======
sampleName = 'ICSV13';
model = 'inputs/tableS2_forward_ICSV13.csv'; % Forward models.
measurements = 'inputs/tableS1_measurements_ICSV13.csv'; % Measurements

% ====== Data type ======
dataFormat = 0; % What type of data do you have? 1 = all measurements (InputA). 0 = mean and std. of variables (InputB).
keep_cols = [1,2,3,4,5,6,7,8,9,10,11]; % Give columns of the model to include. Must include T and P.

% ====== Range of values (only applicable if dataFormat = 0) ======
sd = 0.5; % Range of isopleth values as a multiple of standard deviations from the mean 


%%% CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Functions_NO_EDIT.run_isopleth(sampleName, model, measurements, dataFormat, keep_cols, sd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%