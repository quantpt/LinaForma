function plotResiduals(mu,sigma,model_prediction,obs,dataFormat,variables,sampleName)

% Create distribution of measurements for plotting
if dataFormat == 0
    for i = 1:length(variables)
        x = linspace(mu(i)-3*sigma(i), mu(i) + 3*sigma(i), 1000); 
        observed_dist(i,:) = normpdf(x, mu(i), sigma(i));
        x_axis(i,:) = x;
    end
else
    for i = 1:length(variables)
        tmp = obs(:,i);
        [f,x] = ksdensity(tmp);
        observed_dist(i,:) = f;
        x_axis(i,:) = x;
    end
end

% Figure settings
fig = figure(1);
row = ceil(length(variables)/3);
t = tiledlayout('flow');
set(fig,'Units','centimeters')
set(fig,'Position',[0 0 0.9*21 row*1/5*29])

% Residuals plot for each variable at a time
for i = 1:length(variables)

    % Plot distribution and model prediction
    nexttile
    plot(x_axis(i,:),observed_dist(i,:)); hold on
    plot([model_prediction(i),model_prediction(i)],[0 max(observed_dist(i,:))],'r--','LineWidth',2)

    % Plot observation bounds
    if dataFormat == 0
        plot([(mu(:,i)-2*sigma(:,i)) (mu(:,i)-2*sigma(:,i))],[0 max(observed_dist(i,:))])
        plot([(mu(:,i)+2*sigma(:,i)) (mu(:,i)+2*sigma(:,i))],[0 max(observed_dist(i,:))])
    else
        plot([min(obs(:,i)) min(obs(:,i))],[0 max(observed_dist(i,:))])
        plot([max(obs(:,i)) max(obs(:,i))],[0 max(observed_dist(i,:))])        
    end

    % Give plot a title
    mnT = round(mu(i),3,'significant');
    sigT = round(2*sigma(i),2,'significant');
    t = append('Mean = ',string(mnT),' ± ',string(sigT),' (2σ)');
    title(t);
    ylabel('P.D.E.')
    xlabel(string(variables(i)))
end

% Save figure
name = "FIGURES/L1inv_residuals_" + sampleName + ".svg";
saveas(fig,name);


end