#!/bin/bash

#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --job-name=ChIPseqProcessing
#SBATCH --output=ChIPseqProcessing.output
#SBATCH --error=ChIPseqProcessing.error

module load Miniconda3/4.7.10
source activate jbrowse2
module load R/3.6.0-foss-2018b
module load SAMtools/1.3.1-foss-2016b

DIR="/homes/users/mcoronado/scratch/temp/chip-seq_atac-seq"

while read file
	do
		bigWig=$(echo "$file" | cut -f4 | sed "s/.bigwig//g" )
		bigWigID=$(echo $bigWig | sed 's/_/./1' | cut -f1 -d '.')

		awk -v OFS='\t' '$3 = $3 OFS "."' ${DIR}/ChIP-seq_bed/${bigWigID}.bed > ${DIR}/ChIP-seq_bed/${bigWigID}.input.bed
		assembly=$(echo "$file" | cut -f1 )

		bedmap --echo --mean ${DIR}/chromsizes/${assembly}.chrom.window.bed ${DIR}/ChIP-seq_bed/${bigWigID}.input.bed | tr '|' '\t' | awk -v OFS='\t' '$3 = $3 OFS "."' | sed 's/NAN/0.000000/' | cut -f1,2,3,5  > ${DIR}/ChIP-seq_bed/${bigWigID}.mean.bed
		
		bedGraphToBigWig ${DIR}/ChIP-seq_bed/${bigWigID}.mean.bed ${DIR}/chromsizes/${assembly}.chrom.sizes ${DIR}/ChIP-seq_bw/${bigWigID}.bw
done < <(grep bigwig ${DIR}/data.tsv | grep ChIP-seq)
