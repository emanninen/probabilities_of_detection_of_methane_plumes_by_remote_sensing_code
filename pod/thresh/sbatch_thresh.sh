#!/bin/bash
#
#SBATCH -n 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-11:11          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH --partition=serial_requeue,shared,huce_cascade 
#SBATCH --mem-per-cpu=64000
#SBATCH --job-name thresh
#SBATCH --mail-user o
#SBATCH --mail-user emanninen@g.harvard.edu
#SBATCH --mail-type BEGIN
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
#SBATCH -o /n/home00/emanninen/phd/emiss_dist/scripts/pod/thresh/out_thresh/out.txt
#
echo slurm job id: ${SLURM_JOB_ID}
echo slurm job name: ${SLURM_JOB_NAME}

singularity exec ~/geospatial_latest.sif Rscript /n/home00/emanninen/phd/emiss_dist/scripts/pod/thresh/thresh.r
