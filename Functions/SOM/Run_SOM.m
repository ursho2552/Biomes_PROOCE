function Run_SOM(neurons,epochs,directory_data,directory_out,file_out,iter)
% Funtion to run SOM given the number of neurons, and epochs



%{
Parameters:
    neurons (int): Number of neurons to use, where nuerons^2 represents the
       total number of neurons
    epoch (int): Number of iterations to train the SOM
    directory_data (str) : Path to file where observations are stored
    directory_out (str): Path to folder to store trained SOM and error
    file_out (str): Name of the output file
    iter (int): Iterator in the for loop calling this function
 
 Output
     None. This function saves the trained SOM and errors to the specified
     file
%}

% =========================================================================
% Load data
% e.g. the one named "No_nan_phyto_simple"
% =========================================================================

load(directory_data)



test_dim = [neurons neurons];
test_epoch = epochs;

tic
[classes, net] = My_SOM( No_nan_phyto_simple, test_dim(1),...
    test_dim(2),test_epoch,'mandist' );

cd(directory_out)
%save each SOM run
save(horzcat(file_out,'_',num2str(iter)),'classes','net','iter')
computing_time(i) = toc

tic
[ qe, te, total_error ] = get_total_error( No_nan_phyto_simple,classes, net,'mandist' ); 
%store qe and te to later compare their fraction to the total error
toc

%save the error
save(horzcat(file_out,'_error_',num2str(iter)),'total_error','qe','te','iter','computing_time')
 

exit
end