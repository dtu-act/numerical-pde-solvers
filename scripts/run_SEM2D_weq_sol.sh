#!/bin/bash

#SBATCH -t 05:00:00
#SBATCH --mem=64gb
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=32

#SBATCH -o /users/nborrelj/data/nborrelj/logs/matlab_weq_data%j.out
#SBATCH -e /users/nborrelj/data/nborrelj/logs/matlab_weq_data%j.err
#SBATCH --job-name=matlab_weq_data
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nikolas_borrel-jensen@brown.edu

module load gcc/10.2
module load cuda/11.1.1
module load cudnn/8.2.0
module load python/3.9.0
module load matlab/R2021a

matlab-threaded -nodisplay -r "main_weq_GRF_SEM2D ; exit"