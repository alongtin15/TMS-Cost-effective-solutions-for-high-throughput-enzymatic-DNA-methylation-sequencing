#################################
#SHELL COMMANDS FOR MAKING BSSEQ#
#################################

#create a file for these commands (ie. bsseq.slurm)

#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00

module load GCC/10.2.0
module load OpenMPI/4.0.5
module load R/4.0.5

Rscript making_bsseq.R

#write out and exit the nano
sbatch --mem=100GB --mail-user=email@address --mail-type=ALL bsseq.slurm