%% Test error depending on definition of centroids

%==========================================================================
% In this script we test different approaches to calculate the value that
% is most representative for each cluster
%==========================================================================


cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/01NeuronsError/')
load('Single_run_11.mat')

%load help variables
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/')
load('HelpVariables.mat')
load('Area_map.mat')

%Load simple sort data
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/00Probabilities/')
load('Simple_sort_Data.mat')

%%

tic
original_weights = net.IW{1};
%original_weights = bsxfun(@minus,original_weights,mean(original_weights));
[coeff,score,latent,tsquared,explained] = pca(original_weights);



[r,c] = find(latent > 1)
clearvars latent
%[coeff,score,latent,tsquared,explained]
[coeff,score,~,~,explained] = pca(original_weights,'NumComponents',r(end));
toc



%calculate all, from 2 to 20 biomes
n_clusters = [2:20];
difference_freq = n_clusters.*NaN;
for i = 1:length(n_clusters)

    [raw_monthly_maps, new_weights_freq,new_classes] = Calculate_biomes(net, classes,...
        No_nan_phyto_simple, n_clusters(i),coeff);
  
    
    obs_freq = No_nan_phyto_simple(:,5:end-1) - new_weights_freq(new_classes(:,3),:);
    
    difference_freq(i) = mean(mean(abs(obs_freq),1,'omitnan'),'omitnan');
    

end

%save data
cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/04CentroidDefinition')

if isfile('Centroid_metrics.mat')
    disp('File already exists!')
else
    save('Centroid_metrics','difference_freq')
end



cd(folder_main)
figure
hold on
plot(n_clusters,difference_freq)
plot(n_clusters,difference_single)
grid on
xlabel('Number of clusters')
ylabel('Mean difference between observation and centroid') 

legend('Frequency weighted','Unique weighted')










