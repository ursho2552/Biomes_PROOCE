function Removed_SOM(fr)
% Function traines SOM after removing a fraction fr of the features


    %load original data
    cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/')
    load('Simple_sort_Data.mat')

    %Define optimal setup 
    optimal_dim = [31 31];
    optimal_epoch = 200;

    %get number of features
    N = size(No_nan_phyto_simple(:,5:end-1),2);

    %randomly select fr%
    r_rem = randperm(N,round(fr/100*N));

    tmp_phyto = No_nan_phyto_simple;

    %remove the selected observations
    tmp_phyto(:,r_rem+4) = [];

    %train a SOM
    [removed_classes, removed_net] = My_SOM( tmp_phyto, optimal_dim(1),...
        optimal_dim(2),optimal_epoch,'mandist' );


    str1 = horzcat('Removed_SOM_frac_',int2str(fr));
    cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/06Robustness/FeatureLoss/SOMs/')
    save(str1,'tmp_phyto','removed_classes','removed_net', 'r_rem')
    exit
    
end