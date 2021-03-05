#!/bin/sh
module load matlab/R2017b

## declare an array variable with the values for neurons
declare -a neurons=(5,8,11,14,17,20,23,25,27,29,31,33,35,50,70)
#declare -a epochs=(1, 5, 10, 20, 50, 100, 200, 300, 400, 500, 700, 1000)

# get length of an array to define number of iterations
arraylength=${#neurons[@]}
#arraylength=${#epochs[@]}

#define other input arguments
#neurons=31
epochs=200
directory_data="/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/"
directory_out_neuron="/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/01Neurons_error/"
directory_out_epoch="/net/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/02Epoch_error/"

file_out="Single_run"


# use for loop to read all values and indexes
for (( i=1; i<${arraylength}+1; i++ ));
do
	  bsub -n 36 -W "24:00" matlab -nodisplay -nojvm -r "Run_SOM(${neurons[$i-1]},$epochs,$directory_data,$directory_out_neuron,$file_out,$i)"
      #bsub -n 36 -W "24:00" matlab -nodisplay -nojvm -r "Run_SOM($neurons,${epochs[$i-1]},$directory_data,$directory_out_epoch,$file_out,$i)"

done



