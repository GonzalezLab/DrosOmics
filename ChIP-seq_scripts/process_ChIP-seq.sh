#!/bin/bash

#SBATCH --partition=normal
#SBATCH --cpus-per-task=6
#SBATCH --job-name=ChIPseqProcessing
#SBATCH --output=ChIPseqProcessing.output
#SBATCH --error=ChIPseqProcessing.error


module load Miniconda3/4.7.10
source activate jbrowse2
module load R/3.6.0-foss-2018b
module load SAMtools/1.3.1-foss-2016b

#cp -r /homes/users/mcoronado/gonzalez_lab/browserData/ChIP-seq /homes/users/mcoronado/scratch/temp/chip-seq_atac-seq
#cp -r /homes/users/mcoronado/gonzalez_lab/browserData/ATAC-seq /homes/users/mcoronado/scratch/temp/chip-seq_atac-seq

DIR=/homes/users/mcoronado/scratch/temp/chip-seq_atac-seq
DATA=/homes/users/mcoronado/gonzalez_lab/browserData

mkdir ${DIR}/genomeAssemblies/

for assembly in `less ${DIR}/genomeAssembly.txt`
do
	mkdir ${DIR}/genomeAssemblies/${assembly}
# CLEAN FASTA
genomeAssembly=$(awk -v assembly="$assembly" ' assembly == $1 ' ${DATA}/data.tsv | cut -f 5)
cp ${DATA}/GenomeAssemblies/${genomeAssembly} ${DIR}/genomeAssemblies/${assembly}
seqtk subseq ${DIR}/genomeAssemblies/${assembly}/${genomeAssembly} ${DIR}/chrName.lst > ${DIR}/genomeAssemblies/${assembly}/${assembly}.tmp.fasta
samtools faidx ${DIR}/genomeAssemblies/${assembly}/${assembly}.tmp.fasta $(cat ${DIR}/chrName.lst) > ${DIR}/genomeAssemblies/${assembly}/${assembly}.fasta

rm ${DIR}/genomeAssemblies/${assembly}/${assembly}.tmp.fasta
rm ${DIR}/genomeAssemblies/${assembly}/${assembly}.tmp.fasta.fai
rm ${DIR}/genomeAssemblies/${assembly}/${genomeAssembly}
done

mkdir ${DIR}/chip-seq_atac-seq/chromsizes

while read assembly
do
	echo $assembly
	rm -f ${DIR}/chromsizes/${assembly}.chrom.2.size
	faidx ${DIR}/genomeAssemblies/${assembly}/${assembly}.fasta -i chromsizes > ${DIR}/chromsizes/${assembly}.chrom.sizes

	cat ${DIR}/chromsizes/${assembly}.chrom.sizes | awk -v OFS='\t' '$1 = $1 OFS "0"' >  ${DIR}/chromsizes/${assembly}.chrom.2.sizes
	bedtools makewindows -b ${DIR}/chromsizes/${assembly}.chrom.2.sizes -w 10 > ${DIR}/chromsizes/${assembly}.chrom.window.bed

done < <(cut -f1 ${DIR}/data.tsv | sort -u) 

mkdir $DIR/ATAC-seq_bed
while read file
do
	bigWig=$(echo "$file" |cut -f 4 |sed "s/.bigwig//g" )
	bigWigID=$(echo $bigWig | sed 's/_/./1' | cut -f1 -d '.')
	echo $bigWigID $bigWig
	bigWigToBedGraph ${DIR}/ATAC-seq/${bigWig}.bigwig ${DIR}/ATAC-seq_bed/${bigWigID}.bed
#rm -rf ${wDir}/chipSeq/${assembly}/${tissue}/${histone}/results/call-macs2_signal_track/$bigWig.bigwig
done < <(grep bigwig ${DIR}/data.tsv | grep ATAC-seq) 

while read file
do
	bigWig=$(echo "$file" |cut -f 4 |sed "s/.bigwig//g" )
	bigWigID=$(echo $bigWig | sed 's/_/./1' | cut -f1 -d '.')
	echo $bigWigID $bigWig
	bigWigToBedGraph ${DIR}/ChIP-seq/${bigWig}.bigwig ${DIR}/ChIP-seq_bed/${bigWigID}.bed
#rm -rf ${wDir}/chipSeq/${assembly}/${tissue}/${histone}/results/call-macs2_signal_track/$bigWig.bigwig
done < <(grep bigwig ${DIR}/data.tsv | grep ChIP-seq) 

mkdir ChIP-seq_bw
while read file
do
	bigWig=$(echo "$file" | cut -f4 | sed "s/.bigwig//g" )
	bigWigID=$(echo $bigWig | sed 's/_/./1' | cut -f1 -d '.')

	awk -v OFS='\t' '$3 = $3 OFS "."' ${DIR}/ChIP-seq_bed/${bigWigID}.bed > ${DIR}/ChIP-seq_bed/${bigWigID}.input.bed
	assembly=$(echo "$file" | cut -f1 )

	bedmap --echo --mean ${DIR}/chromsizes/${assembly}.chrom.window.bed ${DIR}/ChIP-seq_bed/${bigWigID}.input.bed | tr '|' '\t' | awk -v OFS='\t' '$3 = $3 OFS "."' | sed 's/NAN/0.000000/' | cut -f1,2,3,5  > ${DIR}/ChIP-seq_bed/${bigWigID}.mean.bed
	bedGraphToBigWig ${DIR}/ChIP-seq_bed/${bigWigID}.mean.bed ${DIR}/chromsizes/${assembly}.chrom.sizes ${DIR}/ChIP-seq_bw/${bigWigID}.bw
done < <(grep bigwig ${DIR}/data.tsv | grep ChIP-seq)

mkdir ATAC-seq_bw
while read file
do
	bigWig=$(echo "$file" | cut -f4 | sed "s/.bigwig//g" )
	bigWigID=$(echo $bigWig | sed 's/_/./1' | cut -f1 -d '.')

	awk -v OFS='\t' '$3 = $3 OFS "."' ${DIR}/ATAC-seq_bed/${bigWigID}.bed > ${DIR}/ATAC-seq_bed/${bigWigID}.input.bed
	assembly=$(echo "$file" | cut -f1 )

	bedmap --echo --mean ${DIR}/chromsizes/${assembly}.chrom.window.bed ${DIR}/ATAC-seq_bed/${bigWigID}.input.bed | tr '|' '\t' | awk -v OFS='\t' '$3 = $3 OFS "."' | sed 's/NAN/0.000000/' | cut -f1,2,3,5  > ${DIR}/ATAC-seq_bed/${bigWigID}.mean.bed
	bedGraphToBigWig ${DIR}/ATAC-seq_bed/${bigWigID}.mean.bed ${DIR}/chromsizes/${assembly}.chrom.sizes ${DIR}/ATAC-seq_bw/${bigWigID}.bw
done < <(grep bigwig ${DIR}/data.tsv | grep ATAC-seq)


while read file
do
	assembly=$(echo "${file}" | cut -f1)
	file=$(echo "${file}" | cut -f4 )
	peakID=$(echo $file | sed 's/_/./1' | cut -f1 -d '.')
	gunzip $DIR/ChIP-seq/${file}
	file=$(echo "${file}" | sed 's/.gz//g')
	sort -k1,1 -k2,2n $DIR/ChIP-seq/${file} > $DIR/ChIP-seq/${file}.sort
# narrow peak: http://genome.ucsc.edu/goldenPath/help/examples/bigNarrowPeak.as
bedToBigBed -type=bed4+6 -as=narrowPeak.as $DIR/ChIP-seq/${file}.sort ${DIR}/chromsizes/${assembly}.chrom.sizes ${DIR}/ChIP-seq_bw/${peakID}.bigBed

#rm -rf ${wDir}/chipSeq/${assembly}/${tissue}/${histone}/results/call-macs2_signal_track/$bigWig.bigwig
done < <(grep Peak.gz ${DIR}/data.tsv | grep ChIP-seq)


while read file
do
	assembly=$(echo "${file}" | cut -f1)
	file=$(echo "${file}" | cut -f4 )
	peakID=$(echo $file | sed 's/_/./1' | cut -f1 -d '.')
	gunzip $DIR/ATAC-seq/${file}
	file=$(echo "${file}" | sed 's/.gz//g')
	sort -k1,1 -k2,2n $DIR/ATAC-seq/${file} > $DIR/ATAC-seq/${file}.sort
# narrow peak: http://genome.ucsc.edu/goldenPath/help/examples/bigNarrowPeak.as
bedToBigBed -type=bed4+6 -as=narrowPeak.as $DIR/ATAC-seq/${file}.sort ${DIR}/chromsizes/${assembly}.chrom.sizes ${DIR}/ATAC-seq_bw/${peakID}.bigBed

#rm -rf ${wDir}/chipSeq/${assembly}/${tissue}/${histone}/results/call-macs2_signal_track/$bigWig.bigwig
done < <(grep Peak.gz ${DIR}/data.tsv | grep ATAC-seq)


while read file; do 
	assembly=$(echo "${file}" | cut -f1)
	type=$(echo "${file}" | cut -f2 | sed -e "s/\b\(.\)/\u\1/g")
	tissue=$(echo "${file}" | cut -f3 | sed -e "s/\b\(.\)/\u\1/g")
	peakID=$(echo "${file}" | cut -f4| sed 's/_/./1' | cut -f1 -d '.')
	if echo "${file}" | grep -q bigwig ; then
		format="bigwig"
	elif echo "${file}" | grep -q Peak ; then
		format="narrowPeak"
	fi
	data=$(echo "${file}" | cut -f5)
	echo -e "$assembly\t$type\t$tissue\t$peakID\t$format\t$data"
done < data.tsv