function plotBoot(t_best,p_best,sampleName)

% Calculate distribution evolution
medT = zeros(1, length(t_best));
q1T = zeros(1, length(t_best));
medP = zeros(1, length(t_best));
q1P = zeros(1, length(t_best));
for i = 1:length(t_best)
    if i == 1
        medT(i) = 0;
        q1T(i) = 0;
        medP(i) = 0;
        q1P(i) = 0;
    else
        preceding_values1 = t_best(1:i-1);
        medT(i) = median(preceding_values1);
        q1T(i) = prctile(preceding_values1,25);
        q3T(i) = prctile(preceding_values1,75);
        preceding_values2 = p_best(1:i-1);
        medP(i) = median(preceding_values2);
        q1P(i) = prctile(preceding_values2,25);
        q3P(i) = prctile(preceding_values2,75);
    end
end
medT = medT(2:end); medP = medP(2:end);
q1T = q1T(2:end); q1P = q1P(2:end);
q3T = q3T(2:end); q3P = q3P(2:end);


% Plot results
figX = figure;
x = 2:length(t_best);
nexttile
errorbar(x,medT,(medT-q1T),(q3T-medT)); xlabel('Bootstrap resamples'); ylabel('Temperature (°C)'); title('Median (+IQR) Temperature'); hold on
plot(x,medT,'r-','LineWidth',2);
legend({'IQR','Median'})
nexttile
errorbar(x,medP,(medP-q1P),(q3P-medP)); xlabel('Bootstrap resamples'); ylabel('Pressure (kbar)'); title('Median (+IQR) Pressure'); hold on
plot(x,medP,'r-','LineWidth',2);
legend({'IQR','Median'})
% Save figure
name = "FIGURES/L1inv_bootstrapDistribution_" + sampleName + ".svg";
saveas(figX,name);

end