function [ bmu2 ] = find_2_bmu( data, class, weights, classes, dist_metric)
%This function find the second best matching unit (neuron) based on
%distance

%{
Parameters:
    data (matrix): Features of one observation
    class (vector): Class label of the observation
    weights (matrix) : Weights of the trained neurons
    classes (vector): Vector containing the class labels for each observations
    dist_metric (str): Metric used to calculate distance between neurons
 
 Output:
    bmu2 (float): label of second best matching unit

%}

    %copy occurrence data
    new_data = data; 

    %find the index in weights that contains the BMU
    ind = find(classes == class); 

    [m, n] = size(weights);
    rest_BMU = ones(m,n+1)*NaN;

    for i = 1:m
    %rest_BMU contains all weights and labels of neurons that are not BMU
    rest_BMU(i,:) = [weights(i,:), classes(i)];   

    end

    %delete BMU, so you only have the non-BMU 
    rest_BMU(ind,:) = [];      

    sum_dist = ones( length(rest_BMU(:,end)),2 )*NaN;
    for i =1:length(rest_BMU(:,end))
        if(strcmp(dist_metric,'dist'))
            dist = new_data - rest_BMU(i,1:end-1);
            dist2 = dist.^2;
            dist3 = sum(dist2,2);
            dist4 = sqrt(dist3);
        elseif(strcmp(dist_metric,'mandist'))
            dist= abs(new_data - rest_BMU(i,1:end-1));
            dist3 = sum(dist,2);
            dist4 = dist3;
        end
        sum_dist(i,1) = dist4;
        sum_dist(i,2) = rest_BMU(i,end);

    end

    %need to check if there's only one BMU, if not, choose one at random
    pre_bmu2 =   sum_dist(find(sum_dist(:,1) == min(sum_dist(:,1))),2);
    if (size(pre_bmu2) ~= [1 1]) 
        fprintf('too many \n')
    end
    bmu2 = datasample(pre_bmu2,1);



   

end

