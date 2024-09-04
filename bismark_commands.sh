#####################
#USING BISMARK SUITE#
#####################

#create a file for the script (ie. bismark.sh)


#!/bin/bash
#
#SBATCH --mem=200GB

module load GCC/8.2.0
module load Intel/2019.1.144
module load Bowtie2/2.3.5.1
module load GCC/11.3.0
module load SAMtools/1.18
module load Bismark/0.24.0

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p /path/to/sampleID`
echo "SampleID: ${sampleID}"

input_path=/path/to/trimmed/files

bismark \
	--genome_folder /path/to/genome/folder/with/BisulfitePrepared/genome/ \
	-1 ${input_path}/${sampleID}_R1_trimmed_paired \
	-2 ${input_path}/${sampleID}_R2_trimmed_paired \
	-o /path/to/output/directory/output_file_name \
	--score_min L,0,-0.6 -R 12 \
	--parallel 8

echo "SampleID: ${sampleID}"

bismark_methylation_extractor \
	-p -o /path/to/output/directory/output_file_name \
	--bedGraph --comprehensive --parallel 24 \
    --genome_folder /path/to/genome/folder/with/BisulfitePrepared/genome/ \
    /path/to/file/produced/from/previous/step/


echo "SampleID: ${sampleID}"

coverage2cytosine \
	--merge_CpG --genome_folder /path/to/genome/folder/with/BisulfitePrepared/genome/ \
	-o /path/to/output/directory/output_file_name \
	/path/to/file/produced/from/previous/step/


#write out and exit the nano
sbatch -t 0-48 --mem=200G --array=1-__ --mail-user=email@address --mail-type=ALL <file name>.sh