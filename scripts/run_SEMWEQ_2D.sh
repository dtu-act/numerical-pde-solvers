#!/bin/bash

#BSUB -W 24:00
#BSUB -q hpc
# -- Number of cores requested -- 
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=16GB]"
#BSUB -J matlab2D_ppw4

### -- Notify me by email when execution begins --
#BSUB -B
### -- Notify me by email when execution ends   --
#BSUB -N
### -- Specify the output and error file. %J is the job-id --
#BSUB -o "/work3/nibor/data/logs/matlab_%J.out"
#BSUB -e "/work3/nibor/data/logs/matlab_%J.err"

module load matlab/R2022a

matlab -nodisplay -batch "SEM_WEQ2D ; exit"