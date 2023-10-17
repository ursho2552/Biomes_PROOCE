function Leaky_SOM(neurons, epochs, fr, directory_data, directory_noise, directory_output)
% Function traines SOM after removing a fraction fr of the observations

%{
Parameters:
    neurons (int): Number of neurons
    epochs (int): Number of epochs
    fr (int): Percentage of features that should be excluded from the
        original dataset prior to training. Beware that this function might 
        seem counterintuitive as it uses an external matrix to determine which
        observations will be excluded.
    directory_data (char): Path to observations
    directory_noise (char): Path to noisy data
    directory_output (char): Path for output
 
 
 Output:
    None. This function directly saves the trained SOMs to the hardcoded
    directory (Should be changed at some point.

%}

    %load original data
    %cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/00Probabilities/')
    load(directory_data)
    
    %cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/06Robustness/SpatialTemporalLoss/')
    %load('Noisy_data_ind.mat')
    load(directory_noise)
    
    %Define optimal setup 
    optimal_dim = [neurons neurons];
    optimal_epoch = epochs;

    %get indeces for fraction
    indeces = sum(~isnan(noisy_data_ind(:,1:fr)),2);
    
    noisy_data = No_nan_phyto_simple;
    noisy_data(No_nan_phyto_simple(indeces == 1,1),:) = [];
    
    [classes_noise, noise_net] = My_SOM( noisy_data, optimal_dim(1),...
    optimal_dim(2),optimal_epoch,'mandist' );

    str1 = horzcat('Leaky_SOM_frac_',int2str(fr));
    
    %cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/06Robustness/SpatialTemporalLoss/SOMs/')
    cd(char(directory_output))
    save(str1,'noisy_data','classes_noise','noise_net','fr','indeces')
    
    exit
    
end