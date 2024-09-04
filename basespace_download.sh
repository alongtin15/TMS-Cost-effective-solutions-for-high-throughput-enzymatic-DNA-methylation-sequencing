####################
#BASESPACE DOWNLOAD#
####################

#website for reference: https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview

#go to target directory
wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs" -O /path/to/directory/bs

#make the file executable
chmod u+x /path/to/directory/bs

#authenticate the file
./bs authenticate

#get list of biosamples you are using for the project
./bs list biosamples

#add the list of biosamples to a file in your directory

#create a file to download the biosamples (ie. download.sh)

#!/bin/bash
#

biosampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p /path/to/directory/biosampleID`
echo "BiosampleID: ${biosampleID}"

./bs biosample download -i ${biosampleID} --extension=<file type> -o /path/to/output/directory/

#write out and exit the nano
sbatch -t 0-10 --mem=50GB --array=1-__ --mail-user=email@address --mail-type=ALL <file name>.sh
