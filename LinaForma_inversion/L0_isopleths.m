%% This script will show overlap between different probe data fields and thermobarometers
clear;clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% INPUTS %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== Data ======
sampleName = 'ICSV17';
model = 'inputs/forward_model.csv'; % Forward models.
measurements = 'inputs/InputA.csv'; % Measurements

% ====== Data type ======
dataFormat = 0; % What type of data do you have? 1 = all measurements (InputA). 0 = mean and std. of variables (InputB).
keep_cols = [1,2,3,4,5,6,7,8,9,10,11]; % Give columns of the model to include. Must include T and P.

% ====== Range of values (only applicable if dataFormat = 0) ======
sd = 0.5; % Range of isopleth values as a multiple of standard deviations from the mean 


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
variablesM = model.Properties.VariableNames;
[X,Y] = meshgrid(T,P);
model_data = data(:,3:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PART 2: Input the observations %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in measurements
obs = readtable(measurements,'ReadRowNames',true,'VariableNamingRule','preserve');
obs = obs(:,keep_cols(3:end)-2);
rows = obs.Properties.RowNames;
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
    observations = [mu+sd*sigma;mu-sd*sigma];
else
    observations = [max(obs,[],2);min(obs,[],1)];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PART 3: Find areas of overlap %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize fields
fields = zeros(length(temperature),length(variables));

% Loop over P-T space for each variable.
for i = 1:length(variables)

    % Get range of values for observation
    obVar = observations(:,i);

    % Loop over P-T space to see if model value lies within range of
    % obserations
    for j = 1:length(temperature)
        modVar = model_data(j,i);
        if modVar >= min(obVar) && modVar <= max(obVar)
            fields(j,i) = 1;
        end
    end
end

% Calculate overlap percentage
percent = sum(fields,2)./length(variables) * 100;
max_percent = max(percent,[],"all");

% Find location of max. overlap field
max_percent_field = zeros(length(temperature),1);
max_percent_field(percent == max_percent) = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PART 4: Plot results %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot overlap percentage
percentGrid = griddata(temperature,pressure,percent,X,Y);
fig1 = figure(1);
set(fig1,'Units','centimeters')
set(fig1,'Position',[0 0 0.9*21 0.9*21])
map = Functions_NO_EDIT.viridis;
pcolor(X,Y,percentGrid); shading flat; colormap(map); c = colorbar; c.Label.String = 'Percentage'; hold on
xlabel('Temperature (°C)')
ylabel('Pressure (bars)')
n = size(fields,2);
sole_field = 1/n *100;
st = append('n compositional fields = ',string(n),'  |   1 field = ', ...
    string(sole_field),'%   |   Max = ',string(max_percent),'% (i.e., ', ...
    string(max(sum(fields,2))),' fields)');
title('Percentage of fields in a given area')
subtitle(st)
axis square

% Plot contours
fig2 = figure(2);
set(fig2,'Units','centimeters')
set(fig2,'Position',[0 0 0.9*21 0.9*21])
map =  Functions_NO_EDIT.linspecer(size(model_data,2));
for i = 1:size(model_data,2)
    grid = griddata(temperature,pressure,model_data(:,i),X,Y);
    [R,e] = contour(X,Y,grid,'Color',map(i,:));
    clabel(R,e,'Color',map(i,:));  
    hold on;
end
leg = variables;
legend(leg)
axis square
xlabel('Temperature (°C)')
ylabel('Pressure (bars)')
title('Contour plot')


% Plot fields of different phases separately
fig3 = figure(3);
set(fig3,'Units','centimeters')
set(fig3,'Position',[0 0 0.9*21 0.9*21])
map =  Functions_NO_EDIT.linspecer(size(model_data,2));
loop = 0;
for i = 1:size(fields,2)
    loop = loop + 1;
    phase = griddata(temperature,pressure,fields(:,i),X,Y);
    phase(phase == 0) = NaN;
    phase = phase + loop;
    p1 = pcolor(X,Y,phase); shading flat; hold on; colormap(map)
    p1.FaceAlpha = 0.4;
end
hold off
legend(variables)
axis square
xlabel('Temperature (°C)')
ylabel('Pressure (bars)')
title('Isopleth fields')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PART 5: Save results %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save table and figures
filename1 = "output_variables\percentageOverlap_" + sampleName + ".mat";
save(filename1,'X','Y','percentGrid');

filename2 = "FIGURES/L0_fig1_" + sampleName + ".svg";
saveas(fig1,filename2);

filename2 = "FIGURES/L0_fig2_" + sampleName + ".svg";
saveas(fig2,filename2);

filename3 = "FIGURES/L0_fig3_" + sampleName + ".svg";
saveas(fig3,filename3);
disp('FINISHED')

%%%%%%%%%%%%%%%%%%%%%
%%%% END OF CODE %%%%
%%%%%%%%%%%%%%%%%%%%%

