#!/bin/bash

#BSUB -W 01:00
#BSUB -q hpc
# -- Number of cores requested -- 
#BSUB -n 32
#BSUB -R "rusage[mem=16GB]"
#BSUB -J matlab_2D

### -- Notify me by email when execution begins --
#BSUB -B
### -- Notify me by email when execution ends   --
#BSUB -N
### -- Specify the output and error file. %J is the job-id --
#BSUB -o "/work3/nibor/1TB/logs/matlab_%J.out"
#BSUB -e "/work3/nibor/1TB/logs/matlab_%J.err"

module load gcc/10.2
module load matlab/R2022b

matlab-threaded -nodisplay -r "SEM_WEQ2D ; exit"