function Removed_SOM(neurons, epochs, fr, directory_data, directory_output)
% Function traines SOM after removing a fraction fr of the features

%{
Parameters:
    neurons (int): Number of neurons
    epochs (int): Number of epochs
    fr (int): Percentage of features that should be excluded from the
        original dataset prior to training
    directory_data (char): Path to observations
    directory_output (char): Path for output
 
 Output:
    None. This function directly saves the trained SOMs to the hardcoded
    directory (Should be changed at some point.

%}

    %load original data
    load(directory_data)

    %Define optimal setup 
    optimal_dim = [neurons neurons];
    optimal_epoch = epochs;

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
    
    %cd('/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/06Robustness/FeatureLoss/SOMs/')
    cd(char(directory_output))
    save(str1,'tmp_phyto','removed_classes','removed_net', 'r_rem')
    exit
    
end