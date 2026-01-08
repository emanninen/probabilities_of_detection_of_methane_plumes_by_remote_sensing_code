#!/bin/bash
#
#SBATCH -n 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 1-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH --partition=serial_requeue,shared,huce_cascade 
#SBATCH --mem-per-cpu=16000
#SBATCH --job-name cr
#SBATCH --mail-user o
#SBATCH --mail-user emanninen@g.harvard.edu
#SBATCH --mail-type BEGIN
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
#SBATCH -o /n/home00/emanninen/phd/emiss_dist/scripts/pod/mair_cont_rel/out_cont_rel.txt  
#
echo 'slurm job id:' $SLURM_JOB_ID
echo 'slurm job name:' $SLURM_JOB_NAME

cont_rel_dir="/n/home00/emanninen/phd/emiss_dist/scripts/pod/mair_cont_rel/"
singularity exec ~/geospatial_latest.sif Rscript ${cont_rel_dir}di_detection.r
