function [ Merger_distance, orig,T] = get_dendrogram(metric,...
    data, labels, numclass,link)
% Function to cluster and plot dendrogram of given observations

%{
Parameters:
    metric (str): Metric used to measure the distance between the different
        componenets. See linkage function:
        (https://ch.mathworks.com/help/stats/linkage.html?searchHighlight=linkage&s_tid=srchtitle)
    data (matrix): Observations to be clustered 
    labels (vector): Labels of the observations
    numclass (int): Number of clusters
    link (str): Linkage method used for clustering. See linkage function:
        (https://ch.mathworks.com/help/stats/linkage.html?searchHighlight=linkage&s_tid=srchtitle)
 
 Output:
    Merger_distance (matrix): Matrix with the distances between clusters
    orig (vector): Vector containing unsorted labels (only used for
        plotting)
    T (matrix): Cluster label of each observation in data after defining
        numclass clusters


%}

%% Select method and perform clustering

    Z = linkage(data,link,metric);
    Merger_distance = Z;
    T = cluster(Z,'maxclust',numclass);
%% Print figures


    figure
    [h,nodes, orig] = dendrogram(Z,0);
    dendro_clusters = gcf;
    figure(dendro_clusters);
    hold on;
    new_labels = labels(orig);
    xticklabels(new_labels);
    hold off;



end



