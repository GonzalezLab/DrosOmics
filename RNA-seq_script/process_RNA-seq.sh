#!/bin/bash

#SBATCH --partition=normal
#SBATCH --cpus-per-task=6
#SBATCH --job-name=bamCoverage
#SBATCH --output=bamCoverage.outputZ
#SBATCH --error=bamCoverage.error

module load Miniconda3/4.7.10
source jbrowse2
module load SAMtools/1.3.1-foss-2016b

DIR="/homes/users/mcoronado/scratch/temp/map"

while read bam
	do
		bamFile=$(echo "$bam" | cut -f4 )
		samtools index ${DIR}/RNA-seq/BAM/${bamFile}.bam
		bamCoverage --bam ${DIR}/RNA-seq/BAM/${bamFile}.bam -o ${DIR}/RNA-seq/RNA-seq/${bamFile}.bw --outFileFormat bigwig --binSize 10 -p 6 --verbose --normalizeUsing BPM
done < data.tsv
# 