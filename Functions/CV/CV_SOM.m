function CV_SOM(fr,Num)
%This traines a SOM for the cross-validation experiment based

%{
Parameters:
    fr (int): Fraction of data leave-out fraction. Note that the value
    is only a name and not the fraction itself. This is not ideal, but will
    leave it as is for the moment
    num (int): ID of fold combination used
 
 Output:
    None. This function saves the output directly to a specified directory.
    Note that the directory is hardcoded

%}

%load data
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities')
load('Simple_sort_Data.mat')

optimal_dim = [31 31];
optimal_epoch = 200;


%load partitions
cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/03CrossValidation/Folds/')
str = horzcat('Data_partitioning_cross_validation_fr_',int2str(fr),'_Seed_7.mat')
load(str)

%select which  partition out of the 10
training_flag = idxtrain(:,Num);


noisy_data = No_nan_phyto_simple(training_flag == 1,:);



[classes_noise, noise_net] = My_SOM( noisy_data, optimal_dim(1),...
    optimal_dim(2),optimal_epoch,'mandist' );

str1 = horzcat('CV_Fr_',int2str(fr),'_Num_',int2str(Num));

cd('/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/03CrossValidation/SOMs')
save(str1,'classes_noise','noise_net','fr','Num','training_flag')
exit

end