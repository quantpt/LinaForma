%%% This script will collate all the text files into a single Excel table
clear;clc;

%%%%%%%%% INPUTS %%%%%%%%%
% Files
domino_folder = 'dom_files';
output_file = 'forward_models.csv';
assemblage_file = 'pixa.csv';

% Pix info
min_T = 650;
max_T = 750;
min_P = 3000;
max_P = 8000;
ix = 100;
iy = 100;


%%%%%%%%%%%%%%%%%%%%% CODE %%%%%%%%%%%%%%%%%%%%
%%%% BEST NOT TO ALTER UNLESS YOU ARE SURE %%%
assemblage_file = append(domino_folder,'/',assemblage_file);
datadir = domino_folder;
fn = dir(datadir + "\*_*");
colNames = zeros(1,length(fn));
precision = ix * iy;
all_phases = zeros([precision,length(fn)+2]);
temp = min_T:(max_T-min_T)/(ix - 1):max_T;
temperature = repmat(temp,iy).';
temperature = temperature(:,1);
pressure = min_P:(max_P-min_P)/(iy - 1):max_P;
pressure = repelem(pressure,ix);
all_phases(:,1) = temperature;
all_phases(:,2) = pressure;
for i = 1:length(fn)
    file = fullfile(fn(i).folder,fn(i).name);
    input_data = readmatrix(file);
    variable = fn(i).name;
    for ii = 1:precision
        n = find(input_data(:,1) == ii);
        if ~isempty(n)
            n_all(ii) = n;
            all_phases(ii,i+2) = input_data(n,2);
        end
    end
end

colNames = struct2cell(fn);


% INPUT ASSEMBLAGE DATA
file = fullfile(fn(i).folder,'assemblage');
input_data = readmatrix(file);
variable = fn(i).name;
col_ass = input_data(:,2);

% Variance calculator
[v,ass] = variance_NO_EDIT(readtable(assemblage_file));
for i = 1:length(col_ass)
    assemblage(i,1) = num2cell(col_ass(i,1));
    value = col_ass(i,1);
    assemblage(i,2) = ass(value,1);
    assemblage(i,3) = num2cell(v(value,1));
end
assemblage = cell2table(assemblage,'VariableNames',["Field","Assemblage","Phases"]);
col = ['Temperature','Pressure',colNames(1,:)];
sTable = array2table(all_phases,'VariableNames',col);
sTable = [sTable(:,[1 2]),assemblage,sTable(:,3:end)];
writetable(sTable,output_file);
disp('FINISHED')