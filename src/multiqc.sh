#!/bin/bash
#SBATCH -p markov_cpu # partition name
#SBATCH -t 0-24:00 # hours:minutes for job
#SBATCH -c 8 # number of cores requested (>= runThredN number)
#SBATCH --mem 64G
#SBATCH --job-name yxz3103_multiqc # Job name
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yxz3103@case.edu
#SBATCH -o multiqc_%j.out # File to which standard out will be written
#SBATCH -e multiqc_%j.err # File to which standard err will be written

if [[ -z "$1" || -z "$2" ]]; then
    echo "Usage: sbatch $0 <ALIGN_DIR> <FASTQC_REPORTS_DIR>"
    exit 1
fi

ALIGN_DIR=$1
FASTQC_REPORTS_DIR=$2

module load Python

pip install local multiqc

~/.local/bin/multiqc "${ALIGN_DIR}"/* "${FASTQC_REPORTS_DIR}"/*
