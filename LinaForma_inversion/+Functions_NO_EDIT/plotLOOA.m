function plotLOOA(dT,dP,variables,tMed,pMed,temperature,pressure)
    figure; hold on
    cmap = jet(length(dP));   % or use parula, jet, hsv, etc.
    markerStyles = {'o','s','d','^','h'};
    nMarkers = numel(markerStyles);
    for fi = 1:length(dP)
        % Cycle through markers
        mkr = markerStyles{mod(fi-1, nMarkers) + 1};
        plot(dT(fi), dP(fi), mkr, ...
            'MarkerSize', 8, ...
            'MarkerFaceColor', cmap(fi,:), ...
            'MarkerEdgeColor', 'k', ...
            'DisplayName', variables{fi});
    end
    % Plot Best-fit
    plot(tMed, pMed*1000, 'k*','LineWidth',1 , 'MarkerSize', 10, 'DisplayName','Best-fit');
    xlabel('Temperature');
    ylabel('Pressure');
    title('Leave-One-Out bestfit P-T');
    grid on;
    legend
    xlim([min(temperature) max(temperature)])
    ylim([min(pressure) max(pressure)])
    axis square
end