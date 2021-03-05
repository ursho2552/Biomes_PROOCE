function [ qe, te, total_error ] = get_total_error( data,...
    classes, net ,dist_metric)
% Number of neurons is calculated analogously to the method of Fendereski 
% et al., 2014. This is done by performing a series of SOM training runs 
% with different number of neurons. The quality is then assessed based on 
% the average quantization and topological errors (Uriarte and Martin 2005 
% and Kohonen 2000).

%{
Parameters:
     data (matrix): Complete data of observations and features
     classes (vector): Vector containing the class labels for each observations
     net (int) : Trained Self-organizing map
     dist_metric (str): Metric used to calculate distance between neurons
 
 Output
     qe (float): Quantization error based on the average norm of the differences
     between the obervation and the neuron weights across all neurons
     te (float): Topological error
     total_error (float): Sum of quantization and topological error
%}



    %net.IW denotes the weights for each input on every neuron
    refvec = net.IW{1}; 

    %creates a matrix with the weights correpsonding to all inputs
    refvecmat = refvec(classes,:); 
                       
    if(strcmp(dist_metric,'dist'))
        sub = data - refvecmat;

        sub = sub.^2;
        sub = sqrt(sum(sub,2));
        sub = [sub,classes];

        n1 = unique(classes);
        qe =ones(n1(end),1)*NaN;
        for i = 1:n1(end)
            sub_tmp = sum(sub(sub(:,2) == i,1));
            qe(i) = sub_tmp/length(sub(sub(:,2) == i));
        end

        qe(isnan(qe)) = [];
        qe = mean(qe);
    
    elseif(strcmp(dist_metric,'mandist'))
        sub = abs(data - refvecmat);
        sub = sum(sub,2);

        sub = [sub,classes];

        n1 = unique(classes);
        qe =ones(n1(end),1)*NaN;
        for i = 1:n1(end)
            sub_tmp = sum(sub(sub(:,2) == i,1));
            qe(i) = sub_tmp/length(sub(sub(:,2) == i)); %or just mean()
        end

        qe(isnan(qe)) = [];
        qe = mean(qe);
    end
    

    %Calculate the second BMU by finding the second smallest distance between
    %observation and neuron

    [n, m] = size(refvec);
    weights_with_label = ones(n,m+1)*NaN;
    for i = 1:size(refvec,1)
        weights_with_label(i,:) = [refvec(i,:),i];
    end

    bmu2 = ones(size(data,1),1)*NaN;%vector containing the second BMU

    for i = 1:size(data,1)
        bmu2(i) = find_2_bmu( data(i,:),classes(i), refvec, 1:length(refvec),dist_metric);
    end


    %Calculate the function u(x) which tells us if the first BMU and second BMU
    %are adjecent (use classes and bmu2)

    weights = refvec;
    [numNeurons, ~] = size(weights);
    neighbors = sparse(tril(net.layers{1}.distances <= 1.001) - eye(numNeurons));
    full_neighbors = full(neighbors);

    n = size(bmu2,1);
    adjacent = ones(size(bmu2,1),1)*NaN;
    for i = 1:n
      if (classes(i) < bmu2(i))
          %then search in full_neighbors using bmu2 as the first index
          adjacent(i) = full_neighbors(bmu2(i),classes(i));
      else
          %use class as the first index
          adjacent(i) = full_neighbors(classes(i),bmu2(i));
      end
    end
    adjacent(adjacent == 0) = 2;
    adjacent(adjacent ==1) = 0;
    adjacent = adjacent/2;
    te = sum(adjacent)/size(adjacent,1);

    %Get total error out of quantization and topological error
    total_error = te + qe;
 
 
end

