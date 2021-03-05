function [ new_classes] = reduce_classes( T, classes)
%takes the T vector form DaviesBouldinDendrogram, which tells us which old
%neurons correspond to a new neuron. It also takes the old classes
%determined by the SOM. This function reassigns the old neuron labels to
%new labels (new_classes), furthermore it gives a matrix where we can see
%which old neruon labels are now summarized by the new ones (new_labels)


%T is a vector which tells me for every neuron what the new label should be
%classes is a matrix (x by 3) which tells me the label of each observation

    new_classes = classes;
    new_classes(:,3) = NaN;

    for i = 1:size(T,1)
        new_classes(classes(:,3) == i,3) = T(i);
    end

end

