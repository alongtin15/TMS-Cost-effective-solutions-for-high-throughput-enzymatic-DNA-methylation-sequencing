##############
#TMS READ ME#
##############

This repository contains code written for "Cost-effective solutions for high-throughput enzymatic DNA methylation sequencing."

Most steps were run using Vanderbilt University's computing cluster (ACCRE) or the ASU Sol super computer. Paths to directories and file output names were generalized to allow for individuals to produce file names that make the most sense to them. 


EXTRACTING DNA METHYLATION AND CREATING BSSEQ OBJECTS
1. downloading files from basespace
	basespace_download.sh

2. trimming with trimmomatic
	trimmomatic.sh

3. bismark suite of commands
	bismark_commands.sh

4. creating bsseq objects
	making_bsseq.R
	making_bsseq.slurm

5. bsseq preprocessing in R to get coverage and methylation files
	bsseq_preprocess.R


DIRECT COMPARISON WITH EPIC ARRAY AND WGBS
1. comparing_EPIC_&_TTMS.R


EXTENSION TO NON-HUMAN PRIMATES
1. downloading files from basespace
	basespace_download.sh

2. trimming with Trim Galore
	trimgalore.sh

3. bismark suite of commands
	bismark_commands.sh

4. creating bsseq objects
	making_bsseq.R
	making_bsseq.slurm

5. bsseq preprocessing in R to get coverage and methylation files
	bsseq_preprocess.R

6. direct comparison with RRBS
	comparing_RRBS_&_TTMS.R

POWER ANALYSES
1. EMseq_power_simulations.R
