#!/bin/bash

#module purge
#module load anaconda/python-3.8.8/2021.05

counter=0
while [ $counter -le 500 ]
    do
    python analysis_scripts/plot_heights.py $counter
    ((counter++))
    done
echo All done
