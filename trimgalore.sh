###################
#USING TRIMMOMATIC#
###################

#create a file for this script (ie. trimming.sh)


#!/bin/bash
#
#SBATCH --mem=50GB

fastq_path=/path/to/FASTQ/files/

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p /path/to/directory/sampleID`

trim_path=/path/to/output/directory/

export PATH=$PATH:/path/to/bin/with/TrimGalore/application/

## trim adaptors
trim_galore --paired -j 8 -o ${trim_path} \
    --basename ${sampleID} --gzip ${fastq_path}/${sampleID}__R1.fastq.gz ${fastq_path}/${sampleID}_R2.fastq.gz