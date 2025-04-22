function resampled_data = nonpara_boot(it,obs)

% Assuming your matrix is named 'data'
data = obs; % Example data for demonstration

% Number of measurements
num_measurements = size(data,1);

% Number of variables
num_variables = size(data,2);

% Number of re-samples
num_resamples = it;

% Initialize the re-sampled data matrix
resampled_data = zeros(num_resamples, num_variables);

% Create the re-samples
for i = 1:num_resamples
    for j = 1:num_variables
        resampled_data(i, j) = data(randi(num_measurements), j);
    end
end

end
