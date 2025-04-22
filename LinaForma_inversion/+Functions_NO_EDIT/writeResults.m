function writeResults(meanT,stdT, meanP, stdP, medianT, iqrT1, iqrT2, medianP, iqrP1,iqrP2, score, T, sampleName)

% Open a file for writing
name = "results_" + sampleName + ".txt"; 
fileID = fopen(name, 'w');

% Helper function to write to both terminal and file
writeBoth = @(fmt, varargin) fprintf(fmt, varargin{:}) & fprintf(fileID, fmt, varargin{:});

% Calculate sensitivity percentage
temp_sens = T{7,:}./meanT * 100;
temp_sens = round(temp_sens,2);
pres_sens = T{8,:}./meanP * 100;
pres_sens = round(pres_sens,2);


% Round values to 1C and 100 bar
meanT = round(meanT,0); stdT = round(stdT,0);
medianT = round(medianT,0); iqrT1 = round(iqrT1,0); iqrT2 = round(iqrT2,0);
meanP = round(meanP,2); stdP = round(stdP,2);
medianP = round(medianP,2); iqrP1 = round(iqrP1,2); iqrP2 = round(iqrP2,2);
score(2) = round(score(2),0); score(3) = round(score(3),2);


% Write the overall section
writeBoth('Overall\n');
writeBoth('-------\n');
writeBoth('Mean = %.3g ± %.3g °C, %.3g ± %.3g kbar (2σ)\n', meanT, stdT, meanP, stdP);
writeBoth('Median = %.3g °C (IQR = %.3g-%.3g °C), %.3g kbar (IQR = %.3g-%.3g kbar)\n', medianT, iqrT1,iqrT2, medianP, iqrP1,iqrP2);
writeBoth('X_total (Median) = %.3g\n', score(1));
writeBoth('SE = %.3g °C (%.3g %%), %.3g kbar (%.3g %%)\n', score(2), score(2) / medianT * 100, score(3), score(3) / medianP * 100);
writeBoth('# of variables = %.3g\n', score(4));
writeBoth('# of fitted variables = %.3g\n', score(5));
writeBoth('Model resolution = %.3g °C, %.3g kbar\n', score(6), score(7));
writeBoth('Bootstrap resamples = %.3g\n',score(8));
writeBoth('\n');

% Write the single variable diagnostics section
writeBoth('Single Variable Diagnostics (Median)\n');
writeBoth('------------------------------------\n');
writeBoth('%-12s', '');
for i = 1:numel(T.Properties.VariableNames)
    writeBoth(' %-12s', T.Properties.VariableNames{i});
end
writeBoth('\n');
for i = 1:height(T)
    writeBoth('%-12s', T.Properties.RowNames{i});
    for j = 1:width(T)
        writeBoth(' %-12s', sprintf('%.4g', T{i, j}));
    end
    writeBoth('\n');
end


% Write diagnostics
writeBoth('\n');
writeBoth('Result Summary\n');
writeBoth('--------------\n');
if score(1) > 1
    writeBoth('**** Warning: Fit is poor. Inversion score above 1. Remove ill-fitting variables.\n');
else
    writeBoth('Fit is acceptable. Inversion score below 1.\n');
end
X = T.Properties.VariableNames;
values = T{6, :};
values = round(values,2);
for i = 1:numel(X)
    if score(1) > 1
        if values(i) > 1
            writeBoth('**** Warning: %s has a score above 1 (Score: %.2f).\n', X{i}, values(i));
        end
    else
        if values(i) > 1
            if temp_sens(i) > 2 || pres_sens(i) > 2
                writeBoth('**** Warning: %s has a score above 1 (Score: %.2f). The result is sensitive to this variable (T: %.2f%%; P: %.2f%%). Please remove.\n', X{i}, values(i),temp_sens(i),pres_sens(i));
            else
                writeBoth('**** %s has a score above 1 (Score: %.2f), but result sensitivity is low (T: %.2f%%; P: %.2f%%). Keep variable. \n', X{i}, values(i),temp_sens(i),pres_sens(i));
            end
        end
    end
end

% Close the file
fclose(fileID);

end
