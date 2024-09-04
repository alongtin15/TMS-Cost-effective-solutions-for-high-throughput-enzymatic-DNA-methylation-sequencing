###################
#USING TRIMMOMATIC#
###################

#create a file for this script (ie. trimming.sh)


#!/bin/bash
#
#SBATCH --mem=50GB

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p /path/to/directory/sampleID`
echo "SampleID: ${sampleID}"

input_path=/path/to/fastq/files/

OutputForwardPaired=/path/to/output/directory/${sampleID}_output_file_name_paired_R1
OutputForwardUnpaired=/path/to/output/directory/${sampleID}_output_file_name_unpaired_R1

OutputReversePaired=/path/to/output/directory/${sampleID}_output_file_name_paired_R2
OutputReverseUnpaired=/path/to/output/directory/${sampleID}_output_file_name_unpaired_R2

threads=8

adapter_file=/path/to/trimmomatic/adapter/file/

java -jar /path/to/trimmomatic/adapter/file/trimmomatic-0.39.jar PE \
        -threads ${threads} \
        ${input_path}/fastq_file_name_R1 \
        ${input_path}/fastq_file_name_R2 \
        ${OutputForwardPaired} ${OutputForwardUnpaired} \
        ${OutputReversePaired} ${OutputReverseUnpaired} \
        ILLUMINACLIP:${adapter_file}:2:30:10:8:true HEADCROP:3 TRAILING:10 MINLEN:25

#write out and exit the nano
sbatch -t 0-4 --mem=4G --array=1-_ --mail-user=email@address --mail-type=ALL <file name>.sh