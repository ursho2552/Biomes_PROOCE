#!/bin/sh

# =================================
# 	INSTRUCTIONS
# ================================
# On Euler run as:
# "Euler_run_CV.sh"
#
# ===============================

module load matlab/R2017b

# ==============================
# Define variables
# =============================

neurons=31
epochs=200

#declare string arguments
directory_data=\"/nfs/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/00Probabilities/Simple_sort_Data.mat\"
directory_output=\"/nfs/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/03CrossValidation/SOMs/\"

# don't use \" as this will be concatenated with other strings
directory_cv_partitions="/nfs/kryo/work/ursho/PhD/Projects/Biomes/Scripts/Biomes_PROOCE/Data/03CrossValidation/Folds/"

# ==============================
# Run 10 fold validation
# =============================
declare -a cvs=( 10 5 3 2 )
declare -a fractions=( 10 20 30 50 )

length=${#cvs[@]}

for (( i=0; i <= length; i++ ))
	do
        cv=${cvs[$i]}
        fr=${fractions[$i]}

    	for (( Num=1; Num <= cv; Num++ ))
    		do
            
            # define partition file
            file_partition=\"${directory_cv_partitions}Data_partitioning_cross_validation_fr_${cv}_Seed_7.mat\"
            
            # run Cross Validation
            sbatch -n 1 --time=24:00:00 --mem-per-cpu=37000 --wrap "matlab -nodisplay -nojvm -r 'CV_SOM($neurons, $epochs, $cv, $Num, $directory_data, $file_partition, $directory_output)'"
            
       		
    		done

	done
