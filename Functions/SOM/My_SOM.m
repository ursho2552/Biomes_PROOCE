function [ classes , net] = My_SOM( data, dimension1, dimension2, ep, dist_metric)
%This function creates, trains and applies a SOM to your data (No_nan_data)
%it returns the associated classes (i.e. the labels of the BMUs) and the
%network itself.

%{
Parameters:
     No_nan_data (matrix): Complete data of observations and features
     dimensionX (int): Number of columns and rows in the neuron lattice
     ep (int) : Number of epochs to train the Self-organizing map
     dist_metric (str): Metric used to calculate distance between neurons
 
 Output
     classes (vector): Vector containing the class labels for each observations
     net (SOM object): Trained Self-organizing map
%}


    %ensure that the data is 2D
    data = squeeze(data);

    %Only use the data without lon, lat, month, ID
    x = data(:,5:end-1);
    %transpose to get the correct structure for the Matlab function
    x =  x';

    %specify all parameters
    dimensions = [dimension1 dimension2];
    net = selforgmap(dimensions);
    net.layers{1}.distanceFcn = dist_metric;
    net.trainParam.epochs = ep;
    % net.performFcn = 'mae';

    % Train the Network
    [net,tr] = train(net,x);

    % Test the Network
    y = net(x);

    %get the neuron number/class
    classes = vec2ind(y); 

    classes = classes';

end

