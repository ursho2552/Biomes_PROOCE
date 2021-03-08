function [ new_classes] = reduce_classes( T, classes)
%takes the T vector form DaviesBouldinDendrogram, which tells us which old
%neurons correspond to a new neuron. It also takes the old classes
%determined by the SOM. This function reassigns the old neuron labels to
%new labels (new_classes), furthermore it gives a matrix where we can see
%which old neruon labels are now summarized by the new ones (new_labels)

%{
Parameters:
    T (matrix): Matrix of the cluster correspondance of each observation at
    a specific clustering level
    classes (vector): Labels of the observations
 
 Output:
    new_classes (vector): Labels after clustering

%}

    new_classes = classes;
    new_classes(:,3) = NaN;

    for i = 1:size(T,1)
        new_classes(classes(:,3) == i,3) = T(i);
    end

end

