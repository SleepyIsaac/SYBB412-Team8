#!/bin/bash
#SBATCH -p markov_cpu # partition name
#SBATCH -t 0-4:00 # hours:minutes for job
#SBATCH -c 8 # number of cores requested (>= runThredN number)
#SBATCH --mem 64G
#SBATCH --job-name yxz3103_trim # Job name
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yxz3103@case.edu
#SBATCH -o %j.out # File to which standard out will be written
#SBATCH -e %j.err # File to which standard err will be written

INPUT_DIR="/scratch/markov2/users/yxz3103/project/fastq"
OUTPUT_DIR="/scratch/markov2/users/yxz3103/project/fastq_trimmed"
REPORT_DIR="/scratch/markov2/users/yxz3103/project/trimming_report"

mkdir -p $OUTPUT_DIR

for file in ${INPUT_DIR}/*_1.fastq.gz
do
    base=$(basename ${file} _1.fastq.gz)

    /home/yxz3103/SYBB412/fastp -i ${INPUT_DIR}/${base}_1.fastq.gz -I ${INPUT_DIR}/${base}_2.fastq.gz \
          -o ${OUTPUT_DIR}/${base}_trimmed_1.fastq.gz -O ${OUTPUT_DIR}/${base}_trimmed_2.fastq.gz \
          --detect_adapter_for_pe \
          --trim_front1 5 --trim_front2 5 \
          --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 \
          --qualified_quality_phred 20 --unqualified_percent_limit 20 \
          --n_base_limit 5 --length_required 30 \
          --dedup \
          --thread 4 \
          --html ${REPORT_DIR}/${base}_fastp_report.html \
          --json ${REPORT_DIR}/${base}_fastp_report.json
done
