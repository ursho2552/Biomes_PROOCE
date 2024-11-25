function CV_SOM(neurons, epochs, fr, Num, directory_data, file_partition, directory_output)
%This traines a SOM for the cross-validation experiment based

%{
Parameters:
    neurons (int): Number of neurons
    epochs (int): Number of epochs
    fr (int): Fraction of data leave-out fraction. Note that the value
    is only a name and not the fraction itself. 
    Num (int): ID of fold combination used
    directory_data (char): Path to observations
    file_partition (char): Path to partition of the data
    directory_output (char): Path for output
 
 Output:
    None. This function saves the output directly to a specified directory.
    Note that the directory is hardcoded

%}

%load data
load(directory_data)

optimal_dim = [neurons neurons];
optimal_epoch = epochs;

%load partitions
load(file_partition)

%select which  partition out of the 10
training_flag = idxtrain(:,Num);


noisy_data = No_nan_phyto_simple(training_flag == 1,:);


[classes_noise, noise_net] = My_SOM( noisy_data, optimal_dim(1),...
    optimal_dim(2),optimal_epoch,'mandist' );



cd(char(directory_output))

str1 = horzcat('CV_Fr_',int2str(fr),'_Num_',int2str(Num));
save(str1,'classes_noise','noise_net','fr','Num','training_flag')

exit

end




