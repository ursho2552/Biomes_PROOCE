function [raw_map_monthly, new_weights,new_classes] = Calculate_biomes(net, classes_orig, data, n_clusters,coeff)
%Function to construct biomes from the class labels, the number of biomes,
%and the fraction of global ocean area that the smallest biomes can have


%{
Parameters:
    net (network): Trained SOM
    classes_orig (vector): Class label of all observations
    data (matrix) : Presence/absence data
    n_clusters (int): Number of clusters to use
    coeff (matrix): PCA transformation matrix
 
 Output:
    raw_map_monthly (matrix): clustered neurons
    new_weights (matrix): values of average neurons after clustering
    new_classes (vector): class label of all observations after clustering

%}


%from the original weights, calculate the average neurons
classes = [data(:,2:3),classes_orig];
% 
% original_weights = net.IW{1};
% original_weights = bsxfun(@minus,original_weights,mean(original_weights));
% [coeff,score,latent,tsquared,explained] = pca(original_weights,'NumComponents',8);
%take only 8 using eigenvalue > 1 criteria
original_weights = net.IW{1};
if(length(coeff) > 1)
    dat_net = original_weights*coeff;
else
    dat_net = original_weights;
end


%cluster trained neurons into n_clusters clusters
[~,T,~] = DaviesBouldinDendrogram('cityblock', dat_net,classes_orig,...
     2,0,'weighted',0);
 
%length(unique(T_matrix_monthly_Bio(:,n_clusters-1)))

%change associated neuron for each observation to the label of the
%respective cluster
[new_classes] = reduce_classes(T(:,n_clusters-1),...
    classes);
%empty_clusters = n_clusters - length(unique(new_classes(:,3)));

old_weights = net.IW{1};%

new_weights = NaN(n_clusters,size(old_weights,2));

for i = 1:n_clusters   
    new_weights(i,:) = mean(old_weights(T(:,n_clusters-1) == i,:),1,'omitnan');
    %new_weights(i,:) = mean(old_weights(classes_orig(new_classes(:,3) == i),:),1,'omitnan');
end


%use new_weights to draw dendrogram showing similarity between clusters
%%
    % =====================================================================
    % =================== now transform with pca ==========================
    % =================== 24 Jun 2019 =====================================
    % =====================================================================
    
    new_weights_T = new_weights*coeff;
    
    
[~,~,~] = DaviesBouldinDendrogram('cityblock', new_weights_T,new_classes(:,end),...
    2,0,'weighted',0 );



%% Monthly biomes Raw
%construct annual and seasonal biomes from unchanged distribution, then
%reduce patchiness (manually if needed!!)

raw_map_monthly = prepare2plot([data(:,2:4),new_classes(:,end)]);

end