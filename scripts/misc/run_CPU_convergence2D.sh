#!/bin/sh
# embedded options to bsub - start with #BSUB
# -- our name ---
#BSUB -J convergence_CPU_2D_main
# -- choose queue --
#BSUB -q hpc
# -- specify that we need 2GB of memory per core/slot -- 
#BSUB -R "rusage[mem=8GB]"
# -- Notify me by email when execution begins --
#BSUB -B
# -- Notify me by email when execution ends   --
#BSUB -N
# -- email address -- 
# please uncomment the following line and put in your e-mail address,
# if you want to receive e-mail notifications on a non-default address
#BSUB -u nibor@elektro.dtu.dk
# -- Output File --
#BSUB -o logs/Output_%J.txt
# -- Error File --
#BSUB -e logs/Error_%J.txt
# -- estimated wall clock time (execution time): hh:mm -- 
#BSUB -W 05:00 
# -- Number of cores requested -- 
#BSUB -n 1 
# -- end of LSF options -- 

matlab -nodisplay -r convergence_CPU_2D_main