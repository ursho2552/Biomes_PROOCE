#!/bin/sh
module load matlab/R2017b


# Spatial-Temporal loss and feature loss


    for fr in {1,5,10,20,30}
    do
        echo $fr
        bsub -n 36 -W "24:00" matlab -nodisplay -nojvm -r "Leaky_SOM($fr)"
        bsub -n 36 -W "24:00" matlab -nodisplay -nojvm -r "Removed_SOM($fr)"
        
    done



