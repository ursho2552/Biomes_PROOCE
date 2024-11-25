function [ bmu2 ] = find_2_bmu( data, class, weights, classes, dist_metric)
%This function find the second best matching unit (neuron) based on
%distance

%{
Parameters:
    data (matrix): Features of one observation
    class (int): Class label of the observation
    weights (matrix) : Weights of the trained neurons
    classes (vector): Vector containing the class labels for each observations
    dist_metric (str): Metric used to calculate distance between neurons
 
 Output:
    bmu2 (float): label of second best matching unit

%}


    classes = classes(:);  % Reshape to column vector

    % Remove the current BMU from the list of neurons (leave only potential second BMUs)
    rest_BMU = [weights, classes];  % Concatenate weights and class labels
    ind = find(classes == class);   % Find the index of the current BMU
    rest_BMU(ind, :) = [];       

    sum_dist = NaN( length(rest_BMU(:,end)),2 );
    for i =1:length(rest_BMU(:,end))
        
        if(strcmp(dist_metric,'dist'))
            dist = sqrt(sum((data - rest_BMU(i, 1:end-1)).^2));
            
        elseif(strcmp(dist_metric,'mandist'))
            dist = sum(abs(data - rest_BMU(i, 1:end-1)));
            
        end
        
        sum_dist(i, :) = [dist, rest_BMU(i, end)];

    end

    % Find the second BMU (smallest distance after excluding the first BMU)
    min_dist_value = min(sum_dist(:, 1));  % Get minimum distance
    candidate_bmu2 = sum_dist(sum_dist(:, 1) == min_dist_value, 2);  % Find the second BMU

    % If multiple candidates exist for the second BMU, randomly choose one
    if numel(candidate_bmu2) > 1
        bmu2 = datasample(candidate_bmu2, 1);
    else
        bmu2 = candidate_bmu2;
    end

end

%{
% Find the index of the BMU (best matching unit)

    classes = classes(:);  % Reshape to column vector

    % Remove the current BMU from the list of neurons (leave only potential second BMUs)
    rest_BMU = [weights, classes];  % Concatenate weights and class labels
    ind = find(classes == class);   % Find the index of the current BMU
    rest_BMU(ind, :) = [];          % Exclude the current BMU

    % Initialize distance array for the second BMU search
    num_neurons = size(rest_BMU, 1);
    sum_dist = NaN(num_neurons, 2);

    % Compute distances for each remaining neuron
    for i = 1:num_neurons
        if strcmp(dist_metric, 'dist')  % Euclidean distance
            dist = sqrt(sum((data - rest_BMU(i, 1:end-1)).^2));
        elseif strcmp(dist_metric, 'mandist')  % Manhattan distance
            dist = sum(abs(data - rest_BMU(i, 1:end-1)));
        else
            error('Unsupported distance metric.');
        end
        sum_dist(i, :) = [dist, rest_BMU(i, end)];  % Store distance and neuron label
    end

    % Find the second BMU (smallest distance after excluding the first BMU)
    min_dist_value = min(sum_dist(:, 1));  % Get minimum distance
    candidate_bmu2 = sum_dist(sum_dist(:, 1) == min_dist_value, 2);  % Find the second BMU

    % If multiple candidates exist for the second BMU, randomly choose one
    if numel(candidate_bmu2) > 1
        bmu2 = datasample(candidate_bmu2, 1);
    else
        bmu2 = candidate_bmu
%}
