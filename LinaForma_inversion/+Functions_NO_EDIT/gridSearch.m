function [t_best,p_best,model_misfit] = gridSearch(it,samples,model_data,temperature,pressure)

% Initialize
loop1 = 0;
t_best  = zeros(it,1); p_best = zeros(it,1);

% Loop over each sample in bootstrap re-sampling
for j = 1:size(samples,1)
loop1 = loop1 + 1; model_misfit = zeros(length(temperature),1);   

% Find the sample values
sample_values = samples(j,:);

% Loop over every grid point
for ii = 1:length(temperature)
    
    % Get forward model values
    model_values = model_data(ii,:);

    % Calculate the misfit between the sample value and the forward model 
    model_misfit(ii,:) = Functions_NO_EDIT.misfit(model_values,sample_values);
    
end

% Find best T and P
indexes = find(model_misfit == min(model_misfit));
if length(indexes) > 1; indexes = indexes(1); end
t_best(loop1,:) = temperature(indexes);
p_best(loop1,:) = pressure(indexes);

end