function [Merger_distance, T_all,labels] = DaviesBouldinDendrogram( method,net,classes, r1,...
    index,link, get_figs)
%function to cluster data "net" and plotting the resulting dendrogram
%{
INPUTS:
method: string specifying the distance metric
(spearman,euclidean,cityblock)
net: SOM (index = 1) or data matrix (index ~= 1) containing the data, observation as rows and attributes
as columns
r1 and r2: number of clusters as a sequence, e.g. from 2 to 100
get_figs: boolean flag for printing dendrogram
link: string specifying the clustering link
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
    Z0 = linkage(data_nan,link,method);
    Z = linkage(full_data,link,method);
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
        orig = NaN;
        EVA1 = NaN;
        EVA2 = NaN;
        EVA3 = NaN;
    end
 
%% Extend results to include empty neurons
T_all = ones(size(data,1),r2).*NaN;
c = 1;

for i = 1:length(available_classes)
    if(available_classes(i) == 0)
        T_all(i,1:length(T(c,:))) = T(c,:);
        c = c + 1;
    elseif(available_classes(i) == 1)
        T_all(i,:) = NaN;
    end
end


