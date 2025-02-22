#!/bin/bash
#SBATCH -p markov_cpu # partition name
#SBATCH -t 0-24:00 # hours:minutes for job
#SBATCH -c 8 # number of cores requested (>= runThredN number)
#SBATCH --mem 64G
#SBATCH --job-name yxz3103_STAR_index # Job name
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yxz3103@case.edu
#SBATCH -o %j.out # File to which standard out will be written
#SBATCH -e %j.err # File to which standard err will be written

if [[ -z "$1" || -z "$2" ]]; then
    echo "Usage: sbatch $0 <FASTQ_DIR> <OUTPUT_DIR>"
    exit 1
fi

FASTQ_DIR=$1
OUTPUT_DIR=$2

cd /scratch/markov2/users/yxz3103/project
mkdir -p "$OUTPUT_DIR"

/home/yxz3103/FastQC/fastqc "${FASTQ_DIR}"/*fastq.gz -o "$OUTPUT_DIR"
