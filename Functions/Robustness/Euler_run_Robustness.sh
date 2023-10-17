#!/bin/sh

# =====================================================================
#                             INSTRUCTIONS
# =====================================================================
# On Euler run as:
# "Euler_run_Robustness.sh keyword"
#
# keyword = spatial to calculate robustness to spatio-temporal data loss
# keyword = features to calculate robustness to feature loss
# =====================================================================

module load matlab/R2017b

# ==============================
# Define variables
# =============================

optimal_neurons=31
optimal_epochs=200

#declare string arguments
directory_data=\"/nfs/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/00Probabilities/Simple_sort_Data.mat\"
directory_noise=\"/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/06Robustness/SpatialTemporalLoss/Noisy_data_ind.mat\"

directory_output_leaky=\"/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/06Robustness/SpatialTemporalLoss/SOMs/\"
directory_output_removed=\"/net/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/06Robustness/FeatureLoss/SOMs/\"

# Test robusteness regarding spatio-temporal or feature loss
declare -a fractions=( 1 5 10 20 30 )

for fr in ${fractions[@]}
do

    if [ $1 == spatial ]
    then
        echo $fr
        sbatch -n 1 --time=24:00:00 --mem-per-cpu=37000 --wrap "matlab -nodisplay -nojvm -r 'Leaky_SOM($neurons, $epochs, $fr, $directory_data, $directory_noise, $directory_output_leaky)'"

    elif [ $1 == features ]
    then
        echo $fr
        sbatch -n 1 --time=24:00:00 --mem-per-cpu=37000 --wrap "matlab -nodisplay -nojvm -r 'Removed_SOM($neurons, $epochs, $fr, $directory_data, $directory_output_removed)'"

    else
        echo "Wrong input. Choose either spatial or features"
    fi

done
