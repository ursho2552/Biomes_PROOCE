%% Test error depending on definition of centroids

cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/01NeuronsError/')
load('Single_run_11.mat')

%load help variables
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/')
load('HelpVariable.mat')
load('Area_map.mat')

%Load simple sort data
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/')
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

% cd('Explained_threshold')

%calculate all, from 2 to 20 biomes
n_clusters = [2:20];
difference_freq = n_clusters.*NaN;
difference_single = difference_freq;
for i = 1:length(n_clusters)

    [raw_monthly_maps, new_weights_freq,new_weights_single,new_classes] = Calculate_biomesV2(net, classes,...
        No_nan_phyto_simple, n_clusters(i),coeff);
  
    
    obs_freq = No_nan_phyto_simple(:,5:end-1) - new_weights_freq(new_classes(:,3),:);
    obs_single = No_nan_phyto_simple(:,5:end-1) - new_weights_single(new_classes(:,3),:);
    
    difference_freq(i) = mean(mean(abs(obs_freq),1,'omitnan'),'omitnan');
    difference_single(i) = mean(mean(abs(obs_single),1,'omitnan'),'omitnan');
    

end

%save data
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/04CentroidDefinition')
save('Centroid_metrics','difference_freq','difference_single')

cd(folder5)
figure
hold on
plot(n_clusters,difference_freq)
plot(n_clusters,difference_single)
grid on
xlabel('Number of clusters')
ylabel('Mean difference between observation and centroid') 

legend('Frequency weighted','unweighted')










