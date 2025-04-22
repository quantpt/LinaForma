function plotInv(t_best,p_best,X,Y,P,T,model_misfit,confidence_level,ix,iy,plot_type,boxplots,T_bins,P_bins,sampleName)


%%%% PART5: Plot misfit surface and perform statistics %%%%
% Plot grid
fig1 = figure;
set(fig1,'Units','centimeters')
set(fig1,'Position',[0 0 0.9*21 0.9*21])
plot(X,Y,'ko');
xlabel('Temperature (°C)')
ylabel('Pressure (kbar)')
minY = min(Y,[],'all'); maxY = max(Y,[],'all');
minX = min(X,[],'all'); maxX = max(X,[],'all');
ylim([(minY-0.05*maxY) (maxY+0.05*maxY)])
xlim([(minX-0.05*maxX) (maxX+0.05*maxX)])
axis square
title('Model grid')


% Plot the misfit surface
fig2 = figure;
t = tiledlayout(3,2);
map = +Functions_NO_EDIT.viridis;
set(fig2,'Units','centimeters')
set(fig2,'Position',[0 0 0.9*21 0.9*25])
nexttile([2 2])
res = reshape(log(model_misfit),[ix iy])';
if plot_type == 1
    contourf(X,Y,res);
else
    pcolor(X,Y,res); shading flat
end
colormap(flipud(map)); hold on
c = colorbar;
c.Label.String = 'Log misfit';
axis square
ylabel('Pressure (kbar)')
xlabel('Temperature (°C)')
ylim([min(P) max(P)])
xlim([min(T) max(T)])

% Add the best-fit solutions and plot error elipse
plot(t_best(:,1),p_best(:,1),'k.','MarkerSize',10);

% Standard and mean calculation
mu_T = mean(t_best); mu_P = mean(p_best);
covariance = cov(t_best,p_best);
std_T = std(t_best); std_P = std(p_best);

% Median and IQR
tMed = median(t_best); pMed = median(p_best);
iqrT1 = prctile(t_best,25); iqrT2 = prctile(t_best,75);
iqrP1 = prctile(p_best,25); iqrP2 = prctile(p_best,75);


% Mode calculation
mode(1,1) = median(t_best);
mode(1,2) = median(p_best);
rP3 = quantile(p_best,0.75); rP1 = quantile(p_best,0.25); 
rT3 = quantile(t_best,0.75); rT1 = quantile(t_best,0.25); 


% Plot error ellipse
[eigen_vectors, eigen_values] = eig(covariance);
major_axis_length = sqrt(eigen_values(1, 1) * chi2inv(confidence_level, 2));
minor_axis_length = sqrt(eigen_values(2, 2) * chi2inv(confidence_level, 2));
rotation_angle = atan2(eigen_vectors(2, 1), eigen_vectors(1, 1));
theta = linspace(0, 2 * pi, 100);
x_ellipse = major_axis_length * cos(theta);
y_ellipse = minor_axis_length * sin(theta);
x_rotated = x_ellipse * cos(rotation_angle) - y_ellipse * sin(rotation_angle);
y_rotated = x_ellipse * sin(rotation_angle) + y_ellipse * cos(rotation_angle);
x_rotated = x_rotated + mu_T;
y_rotated = y_rotated + mu_P;
plot(x_rotated, y_rotated, 'r-', 'LineWidth', 3);
plot(mu_T,mu_P,"pentagram",'MarkerFaceColor','yellow','MarkerEdgeColor','k','MarkerSize',20)
plot(mode(1,1),mode(1,2),"pentagram",'MarkerFaceColor','blue','MarkerEdgeColor','k','MarkerSize',20)
ellipse_name = append(string(confidence_level),' uncertainty ellipse');
legend('Misfit map','Monte Carlo points of minimum misfit',ellipse_name,'Mean P-T point','Median P-T point')
title('Grid-search results')


% Round to 1C and 100bar
med_T = round(tMed,0); rT1 = round(iqrT1,0); rT3 = round(iqrT2,0);
med_P = round(pMed,2); rP1 = round(iqrP1,2); rP3 = round(iqrP2,2);


% Plot boxplots
if boxplots == 1
    nexttile
    boxplot(t_best,'Orientation','horizontal')
    s = append('Md = ',string(med_T),' °C (IQR =  ',string(rT1),'-',string(rT3),')');
    title(s); clc;
    xlabel('Temperature (°C)')
    nexttile
    boxplot(p_best,'Orientation','horizontal')
    s = append('Md = ',string(med_P),' kbar (IQR =  ',string(rP1),'-',string(rP3),')');
    title(s); clc;
    xlabel('Pressure (kbar)')

% Or plot histograms
else
    % T
    nexttile
    histogram(t_best,T_bins);
    s = append('Md = ',string(med_T),' °C (IQR =  ',string(rT1),'-',string(rT3),')');
    title(s)
    xlabel('Temperature (°C)')
    ylabel('Number of solutions')
    
    % P
    nexttile
    histogram(p_best,P_bins);
    s = append('Md = ',string(med_P),' kbar (IQR =  ',string(rP1),'-',string(rP3),')');
    title(s)
    xlabel('Pressure (kbar)')
    ylabel('Number of solutions')

end

% Plot results on the contour plot
filenameTest = "output_variables\percentageOverlap_" + sampleName + ".mat";
if exist(filenameTest, 'file')
load(filenameTest);
fig3 = figure;
set(fig3,'Units','centimeters')
set(fig3,'Position',[0 0 0.9*21 0.9*21])
pcolor(X,Y,percentGrid); colormap(map); shading flat; c = colorbar; hold on
c.Label.String = 'Percentage of observations which overlap';
plot(t_best(:,1),p_best(:,1),'k.','MarkerSize',10);
plot(mu_T,mu_P,"pentagram",'MarkerFaceColor','yellow','MarkerEdgeColor','k','MarkerSize',20)
plot(mode(1,1),mode(1,2),"pentagram",'MarkerFaceColor','blue','MarkerEdgeColor','k','MarkerSize',20)
axis square
xlabel('Temperature (°C)')
ylabel('Pressure (kbar)'); hold off
title('Best-fit solutions and overlapping contours')
legend('Contour plot','best-fit solutions','Mean best-fit solution','Median best-fit solution')
name = "FIGURES/L1inv_overlap_" + sampleName + ".svg";
saveas(fig3,name);
end

% Plot grid of results
fig4 = figure;
set(fig4,'Units','centimeters')
set(fig4,'Position',[0 0 0.9*21 0.9*21])
histogram2(t_best(:,1),p_best(:,1),[T_bins,P_bins],'DisplayStyle','tile','ShowEmptyBins','on'); 
colormap(map)
c = colorbar;
c.Label.String = 'Log number of solutions'; hold on
set(gca,'ColorScale','log')
plot(mu_T,mu_P,"pentagram",'MarkerFaceColor','yellow','MarkerEdgeColor','k','MarkerSize',20)
plot(mode(1,1),mode(1,2),"pentagram",'MarkerFaceColor','blue','MarkerEdgeColor','k','MarkerSize',20)
axis square
xlabel('Temperature (°C)')
ylabel('Pressure (kbar)')
title('Best-fit solution 2D histogram')
legend('2D histogram','Mean best-fit solution','Median best-fit solution');


% Save plots and PT solutions
name = "FIGURES/L1inv_modelGrid_" + sampleName + ".svg";
saveas(fig1,name);

name = "FIGURES/L1inv_results_" + sampleName + ".svg";
saveas(fig2,name);

name = "FIGURES/L1inv_solutionDistribution_" + sampleName + ".svg";
saveas(fig4,name);


end