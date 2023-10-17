#/bin/sh

# =====================================================================
#                             INSTRUCTIONS
# =====================================================================
# On Euler run as:
# "Euler_run_SOM.sh keyword"
#
# keyword = neurons to calculate optimal number of neurons
# keyword = epochs to calculate optimal number of epochs
# =====================================================================

module load matlab/R2017b

## declare array variables to loop over
declare -a neuron_array=(5 8 11 14 17 20 23 25 27 29 31 33 35 50 70)
declare -a epoch_array=(1 5 10 20 50 100 200 300 400 500 700 1000)

#define other input arguments
optimal_neurons=31
optimal_epochs=200

#declare string arguments
directory_data=\"/nfs/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/00Probabilities/Simple_sort_Data.mat\"
directory_out_neuron=\"/nfs/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/01Neurons_error/\"
directory_out_epoch=\"/nfs/kryo/work/ursho/Damiano_Presence_data/presence_absence_tables_ensemble_averages/Group_specific_background_approach/Data/02Epoch_error/\"
file_out=\"Single_run\"

COUNTER=1
if [ $1 == neurons ]
then
      #loop for testing optimal number of neurons
      for value in "${neuron_array[@]}"
      do
            sbatch -n 1 --time=24:00:00 --mem-per-cpu=37000 --wrap "matlab -nodisplay -nojvm -r 'Run_SOM($value,$optimal_epochs,$directory_data,$directory_out_neuron,$file_out,$COUNTER)'"
            COUNTER=$[COUNTER + 1]
		break
      done

elif [ $1 == epochs ]
then
      #loop for testing optimal number of epochs
      for value in "${epoch_array[@]}"
      do
            sbatch -n 1 --time=24:00:00 --mem-per-cpu=37000 --wrap "matlab -nodisplay -nojvm -r 'Run_SOM($optimal_neurons,$value,$directory_data,$directory_out_epoch,$file_out,$COUNTER)'"
            COUNTER=$[COUNTER + 1]
      done

else
      echo "Wrong input. Choose either neurons or epochs"
fi





