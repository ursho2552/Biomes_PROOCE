function [Merger_distance, T_all,labels] = DaviesBouldinDendrogram( metric,net,classes, r1,...
    index,link, get_figs)
% Function to cluster the trained neurons and plotting the resulting
% dendrogram

%{
Parameters:
    metric (str): Metric used to measure the distance between the different
        componenets. See linkage function:
        (https://ch.mathworks.com/help/stats/linkage.html?searchHighlight=linkage&s_tid=srchtitle)
    net (Network or matrix): Trained SOM 
    classes (vector): Labels of the observations
    r1 (int): Minimum number of clusters
    index (bool): Flag specifying whether the data is a Matlab network or a
        matrix of observations and features. 1 for network, 0 for matrix
    link (str): Linkage method used for clustering. See linkage function:
        (https://ch.mathworks.com/help/stats/linkage.html?searchHighlight=linkage&s_tid=srchtitle)
    get_figs (bool): Flag to indicate whether or not to print the
        dendrogram. 1 to print dendrogram
 
 Output:
    Merger_distance (matrix): Matrix with the distances between clusters
    T_all (matrix): Cluster label of each observation in net after defining r1
    up to the meximum number of possible clusters
    labels (vector): All possible labels given the size of net

%}


%% Prepare data depending on format
    if(index == 1)
        %get the weights from the network
        data = net.IW{1};
    else
        %take the data as is
        data = net;
    end
% get labels from data matrix, as row index     
    labels = 1:size(data,1);
    r2 = labels(end);
% only cluster neurons that were used and not all of them    
    available_classes = ~ismember(labels,unique(classes));
    full_data = data;
    data_nan = data;
    full_labels = labels;
    data_nan(available_classes ==1,:) = NaN;
    full_data(available_classes== 1,:) = [];
    full_labels(available_classes== 1) = [];
    r2_full = 100;%r2-sum(available_classes);
% cluster data
    Z0 = linkage(data_nan,link,metric);
    Z = linkage(full_data,link,metric);
    Merger_distance = Z;
    T = cluster(Z,'maxclust',r1:r2_full);
%% Print figures
    if(get_figs == 1)
        %
    %EVA1 = evalclusters(full_data,T,'silhouette'); %evaluates optimal number of clusters
    %{
    figure
    hold on
    plot(EVA.InspectedK,EVA.CriterionValues,'k.-')
    title('Silhouette')
    hold off;
    %}
    %EVA2 = evalclusters(full_data,T,'CalinskiHarabasz'); %evaluates optimal number of clusters
    %{
    figure
    hold on
    plot(EVA.InspectedK,EVA.CriterionValues,'k.-')
    title('CalinskiHarabasz')
    hold off;
    %}
    %EVA3 = evalclusters(full_data,T,'DaviesBouldin'); %evaluates optimal number of clusters
    %{
    figure
    hold on
    plot(EVA.InspectedK,EVA.CriterionValues,'k.-')
    title('DaviesBouldin')
    hold off;
    %}
%     %}
%         EVA1 = NaN;
%         EVA2 = NaN;
%         EVA3 = NaN;
        figure
        [h,nodes, orig] = dendrogram(Z,0);
        dendro_clusters = gcf;
        figure(dendro_clusters);
        hold on;
        new_labels = full_labels(orig);
        xticklabels(new_labels);
        hold off;

    else
        %this part is deprecated and does not contribute anything to the
        %current function
        orig = NaN;
        EVA1 = NaN;
        EVA2 = NaN;
        EVA3 = NaN;
    end
 
%% Extend results to include empty neurons
T_all = NaN(size(data,1),r2);
c = 1;

for i = 1:length(available_classes)
    if(available_classes(i) == 0)
        T_all(i,1:length(T(c,:))) = T(c,:);
        c = c + 1;
    elseif(available_classes(i) == 1)
        T_all(i,:) = NaN;
    end
end


