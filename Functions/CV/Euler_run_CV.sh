#!/bin/sh
module load matlab/R2017b





# for 10-fold


    for Num in {1..10}
    do
        echo $Num
        bsub -n 36 -W "24:00" matlab -nodisplay -nojvm -r "CV_SOM(10,$Num)"

    done



# for 5-fold


    for Num in {1..5}
    do
        echo $Num
        bsub -n 36 -W "24:00" matlab -nodisplay -nojvm -r "CV_SOM(20,$Num)"

    done



# for 3-fold

    for Num in {1..3}
    do
        echo $Num
        bsub -n 36 -W "24:00" matlab -nodisplay -nojvm -r "CV_SOM(30,$Num)"

    done



# for 2-fold

    for Num in {1..2}
    do
        echo $Num
        bsub -n 36 -W "24:00" matlab -nodisplay -nojvm -r "CV_SOM(50,$Num)"

    done


