function [ classes ] = Apply_SOM( data, net )
%This function applies a trained some to a dataset

    data = squeeze(data);
    x = data(:,5:end-1);
    x =  x';

    % Apply SOM to data
    y = net(x);

    %get the neuron number/class
    classes = vec2ind(y); 

    classes = classes';
end