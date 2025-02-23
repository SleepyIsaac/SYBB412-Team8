#!/bin/bash
#SBATCH -p markov_cpu # partition name
#SBATCH -t 0-4:00 # hours:minutes for job
#SBATCH -c 8 # number of cores requested (>= runThredN number)
#SBATCH --mem 64G
#SBATCH --job-name yxz3103_feature_counts # Job name
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yxz3103@case.edu
#SBATCH -o fea_%j.out # File to which standard out will be written
#SBATCH -e fea_%j.err # File to which standard err will be writtenfeatureCounts -T 8 -s 2 

if [[ -z "$1" || -z "$2" ]]; then
    echo "Usage: qsub $0 <BAM_DIR> <OUTPUT_DIR>"
    exit 1
fi

BAM_DIR=$1
OUTPUT_DIR=$2

cd /scratch/markov2/users/yxz3103/project

mkdir -p $OUTPUT_DIR

featureCounts -T 8 -s 2 -p -a annotation/Homo_sapiens.GRCh38.113.gtf \
	-o "${OUTPUT_DIR}"/gene_id.counts \
	"${BAM_DIR}"/*.bam
