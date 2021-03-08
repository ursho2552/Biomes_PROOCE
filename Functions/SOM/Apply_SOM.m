function [ classes ] = Apply_SOM( data, net )
%This function applies a trained some to a dataset
%{
Parameters:
    data (vector): vectorized dataset with observations x features
    net (Network or matrix): Trained SOM 
 
 Output:
    classes (vector): Vector with the classes linked to the observations in
        data

%}

    data = squeeze(data);
    x = data(:,5:end-1);
    x =  x';

    % Apply SOM to data
    y = net(x);

    %get the neuron number/class
    classes = vec2ind(y); 

    classes = classes';
end