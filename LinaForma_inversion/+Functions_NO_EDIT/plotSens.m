function plotSens(T_min,T_max,P_min,P_max,variables,T_best,P_best,sampleName)

% Plot T variability
fig1 = figure;
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
T_best = round(T_best,0);
barh(sorted_impact,'stacked'); % Blue bars for positive impact
hold on;
set(gca, 'YTick', 1:length(variables), 'YTickLabel', fliplr(sortedVariables),'Fontsize',12); % Set y-axis labels
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

[~, sortOrder] = sort(total_impact, 'descend');
sortedVariables = variables(sortOrder);
sorted_impact = impact(flipud(sortOrder),:);

% Round to 1C and 100 bar
P_best = round(P_best,2);


% Create tornado plot
barh(sorted_impact,'stacked');
set(gca, 'YTick', 1:length(variables), 'YTickLabel', fliplr(sortedVariables),'Fontsize',12); % Set y-axis labels
xlabel('Pressure (bar)'); % Label x-axis
t = append('Variation relative to ',string(P_best),' bar');
title(t)
lg = legend('-','+','Location','southeast');
xlim([min(min(sorted_impact))-0.1 max(max(sorted_impact))+0.1])

% Save figure
name = "FIGURES/L1inv_sensitivity_" + sampleName + ".svg";
saveas(fig1,name);

end
