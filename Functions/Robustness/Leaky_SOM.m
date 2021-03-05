function Leaky_SOM(fr)
% Function traines SOM after removing a fraction fr of the observations


    %load original data
    cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/')
    load('Simple_sort_Data.mat')
    
    cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/06Robustness/SpatialTemporalLoss/')
    load('Noisy_data_ind.mat')
    
    %Define optimal setup 
    optimal_dim = [31 31];
    optimal_epoch = 200;

    
    %get indeces for fraction
    indeces = sum(~isnan(noisy_data_ind(:,1:fr)),2);
    
    noisy_data = No_nan_phyto_simple;
    noisy_data(No_nan_phyto_simple(indeces == 1,1),:) = [];
    
    [classes_noise, noise_net] = My_SOM( noisy_data, optimal_dim(1),...
    optimal_dim(2),optimal_epoch,'mandist' );

    str1 = horzcat('Leaky_SOM_frac_',int2str(fr));
    nn = fr;
    
    cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/06Robustness/SpatialTemporalLoss/SOMs/')
    save(str1,'noisy_data','classes_noise','noise_net','nn','indeces')
    
    exit
    
end