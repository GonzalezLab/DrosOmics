#!/bin/bash

jBrowseDir="/var/www/html/drosOMICS"
dataDir="/home/marta/Desktop/BROWSER"

tracks="genes TE"

mkdir -p  ${jBrowseDir}/DATA/GENOME/
mkdir -p  ${jBrowseDir}/DATA/GENES/
mkdir -p  ${jBrowseDir}/DATA/TE/
mkdir -p  ${jBrowseDir}/DATA/RNA-seq/
mkdir -p  ${jBrowseDir}/DATA/ChIP-seq/
mkdir -p  ${jBrowseDir}/DATA/ATAC-seq/

cd ${jBrowseDir}

# Create genome, gene and TE tracks
while read assembly
	do
		assemblyCode=$(echo "$assembly" | cut -f1)
		description=$(echo "$assembly" | cut -f3)
		assemblyName=$(echo "$assembly" | cut -f4)
		genomeAssembly=$(echo "$assembly" | cut -f5)
		geneAnnotation=$(echo "$assembly" | cut -f6)
		geneAnnotationFolder=$(echo "$geneAnnotation" | cut -f1 -d'/')
		geneAnnotation=$(echo "$geneAnnotation" | cut -f2 -d'/')
		TEannotation=$(echo "$assembly" | cut -f7)

		#### GENOME
		# Extract only chr 2L, 2R, 3L, 3R and X
		seqtk subseq ${dataDir}/GenomeAssemblies/${genomeAssembly} ${dataDir}/chrName.lst > ${dataDir}/GenomeAssemblies/${assemblyCode}.fasta
		cp ${dataDir}/GenomeAssemblies/${assemblyCode}.fasta ${jBrowseDir}/DATA/GENOME/${assemblyCode}.tmp.fasta

		# Order by chr list
		samtools faidx ${jBrowseDir}/DATA/GENOME/${assemblyCode}.tmp.fasta $(cat ${dataDir}/chrName.lst) > ${jBrowseDir}/DATA/GENOME/${assemblyCode}.fasta

		rm ${jBrowseDir}/DATA/GENOME/${assemblyCode}.tmp.fasta

		bgzip -f ${jBrowseDir}/DATA/GENOME/${assemblyCode}.fasta

		samtools faidx ${jBrowseDir}/DATA/GENOME/${assemblyCode}.fasta.gz

		# Upload genome assembly
		jbrowse add-assembly DATA/GENOME/${assemblyCode}.fasta.gz --load inPlace --name "${assemblyCode}" --displayName "${description}" --type bgzipFasta --force

		#### GENE TRACK
		# Clean liftoff anotation and transform gtf to gff3
		cp ${dataDir}/GeneAnnotations/${geneAnnotationFolder}/${geneAnnotation} ${dataDir}/GeneAnnotations/${geneAnnotationFolder}/${assemblyCode}.gtf

		sed 's/coverage \"\S*\"; //g' ${dataDir}/GeneAnnotations/${geneAnnotationFolder}/${assemblyCode}.gtf | sed 's/sequence_ID \"\S*\"; //g'  | sed 's/extra_copy_number \"\S*\"; //g' | sed 's/copy_num_ID \"\S*\"; //g' | sed 's/partial_mapping \"\S*\"; //g'  | sed 's/low_identity \"\S*\"; //g' > ${dataDir}/GeneAnnotations/${geneAnnotationFolder}/${assemblyCode}.clean.gtf

		gffread ${dataDir}/GeneAnnotations/${geneAnnotationFolder}/${assemblyCode}.clean.gtf -O -F -o ${dataDir}/GeneAnnotations/${geneAnnotationFolder}/${assemblyCode}.clean.gff

		gff3sort.pl ${dataDir}/GeneAnnotations/${geneAnnotationFolder}/${assemblyCode}.clean.gff > ${dataDir}/GeneAnnotations/${geneAnnotationFolder}/${assemblyCode}.tmp.gff3

		perl -p -e "s/;Parent=FBgn[0-9]{7}\n/\n/g" ${dataDir}/GeneAnnotations/${geneAnnotationFolder}/${assemblyCode}.tmp.gff3 >  ${dataDir}/GeneAnnotations/${geneAnnotationFolder}/${assemblyCode}.gff3 

		bgzip -f ${dataDir}/GeneAnnotations/${geneAnnotationFolder}/${assemblyCode}.gff3
		
		tabix -f ${dataDir}/GeneAnnotations/${geneAnnotationFolder}/${assemblyCode}.gff3.gz

		cp ${dataDir}/GeneAnnotations/${geneAnnotationFolder}/${assemblyCode}.gff3.gz ${jBrowseDir}/DATA/GENES/
		
		cp ${dataDir}/GeneAnnotations/${geneAnnotationFolder}/${assemblyCode}.gff3.gz.tbi ${jBrowseDir}/DATA/GENES/

		sed "s/{assemblyName}/${assemblyCode}/g" ${dataDir}/genes_template.json > ${dataDir}/genes_${assemblyCode}.json

		jbrowse add-track-json ${dataDir}/genes_${assemblyCode}.json --update

		#### TE track
		# Transform TE annotation from bed to gff3
		awk -F'\t' '{print $1"\tGonzalezLab\tTE\t"$2"\t"$3"\t.\t"$6"\t.\tTE_id="$4"; Family="$7"; " }' ${dataDir}/TE_annotations_common_name/${TEannotation} > ${dataDir}/TE_annotations_common_name/${assemblyCode}.gff
		
		gff3sort.pl ${dataDir}/TE_annotations_common_name/${assemblyCode}.gff > /${dataDir}/TE_annotations_common_name/${assemblyCode}.TE.gff3
		
		bgzip -f ${dataDir}/TE_annotations_common_name/${assemblyCode}.TE.gff3
		
		tabix -f ${dataDir}/TE_annotations_common_name/${assemblyCode}.TE.gff3.gz
		
		cp ${dataDir}/TE_annotations_common_name/${assemblyCode}.TE.gff3.gz ${jBrowseDir}/DATA/TE/
		
		cp ${dataDir}/TE_annotations_common_name/${assemblyCode}.TE.gff3.gz.tbi ${jBrowseDir}/DATA/TE/

		sed "s/{assemblyName}/${assemblyCode}/g" ${dataDir}/TE_template.json > ${dataDir}/TE_${assemblyCode}.json

		jbrowse add-track-json ${dataDir}/TE_${assemblyCode}.json  --update

done < <(tail -n +2 ${dataDir}/data.tsv)

# Index gene and TE names
jbrowse text-index --attributes=gene_id,geneID,gene_symbol,TE_id --exclude=mRNA,CDS,exon --force

#### RNA-seq track
while read sample
	do
		assemblyCode=$(echo "$sample" | cut -f1)
		condition=$(echo "$sample" | cut -f2)
		tissue=$(echo "$sample" | cut -f3)
		sampleID=$(echo "$sample" | cut -f4)
		cp ${dataDir}/RNA-seq/${sampleID}.bw ${jBrowseDir}/DATA/RNA-seq
		if [ $condition == "-" ] # Normal RNA-seq (no stress)
			then
			sed "s/{sample}/${sampleID}/g" ${dataDir}/RNA_seq_template.json > ${dataDir}/${sampleID}_RNA_seq_template.json
			sed -i "s/{tissue}/${tissue}/g" ${dataDir}/${sampleID}_RNA_seq_template.json
			sed -i "s/{assemblyName}/${assemblyCode}/g" ${dataDir}/${sampleID}_RNA_seq_template.json
			jbrowse add-track-json ${dataDir}/${sampleID}_RNA_seq_template.json --update
		elif [ $condition == "Control" ] # RNA-seq in control conditions
			then
			sed "s/{sample}/${sampleID}/g" ${dataDir}/RNA_seq_control_template.json > ${dataDir}/${sampleID}_RNA_seq_control_template.json
			sed -i "s/{tissue}/${tissue}/g" ${dataDir}/${sampleID}_RNA_seq_control_template.json
			sed -i "s/{assemblyName}/${assemblyCode}/g" ${dataDir}/${sampleID}_RNA_seq_control_template.json
			jbrowse add-track-json ${dataDir}/${sampleID}_RNA_seq_control_template.json --update
		else # RNA-seq under a treatment
			sed "s/{sample}/${sampleID}/g" ${dataDir}/RNA_seq_treatment_template.json > ${dataDir}/${sampleID}_RNA_seq_treatment_template.json
			sed -i "s/{tissue}/${tissue}/g" ${dataDir}/${sampleID}_RNA_seq_treatment_template.json
			sed -i "s/{assemblyName}/${assemblyCode}/g" ${dataDir}/${sampleID}_RNA_seq_treatment_template.json
			sed -i "s/{condition}/${condition}/g" ${dataDir}/${sampleID}_RNA_seq_treatment_template.json
			jbrowse add-track-json ${dataDir}/${sampleID}_RNA_seq_treatment_template.json --update
		fi
done < ${dataDir}/data_RNA-seq.tsv

#### ChIP-seq and ATAC-seq track
while read sample
	do
		assemblyCode=$(echo "$sample" | cut -f1)
		condition=$(echo "$sample" | cut -f2)
		histone=$(echo "$sample" | cut -f3)
		tissue=$(echo "$sample" | cut -f4)
		sampleID=$(echo "$sample" | cut -f5)
		type=$(echo "$sample" | cut -f6)
		experiment=$(echo "$sample" | cut -f7)

		if [ $condition == "-" -a $type == "bigwig" -a $experiment == "ChIP-seq" ] 
			then
			cp ${dataDir}/ChIP-seq/${sampleID}.bw ${jBrowseDir}/DATA/ChIP-seq
			sed "s/{sample}/${sampleID}/g" ${dataDir}/ChIP_seq_template.json > ${dataDir}/${sampleID}_ChIP_seq_template.json
			sed -i "s/{tissue}/${tissue}/g" ${dataDir}/${sampleID}_ChIP_seq_template.json
			sed -i "s/{assemblyName}/${assemblyCode}/g" ${dataDir}/${sampleID}_ChIP_seq_template.json
			sed -i "s/{histone}/${histone}/g" ${dataDir}/${sampleID}_ChIP_seq_template.json
			jbrowse add-track-json  ${dataDir}/${sampleID}_ChIP_seq_template.json --update

		elif [ $condition == "-" -a $type == "narrowPeak" -a $experiment == "ChIP-seq" ] 
			then
			cp ${dataDir}/ChIP-seq/${sampleID}.bigBed ${jBrowseDir}/DATA/ChIP-seq
			sed "s/{sample}/${sampleID}/g" ${dataDir}/ChIP_seq_template_peak.json > ${dataDir}/${sampleID}_ChIP_seq_template_peak.json
			sed -i "s/{tissue}/${tissue}/g" ${dataDir}/${sampleID}_ChIP_seq_template_peak.json
			sed -i "s/{assemblyName}/${assemblyCode}/g" ${dataDir}/${sampleID}_ChIP_seq_template_peak.json
			sed -i "s/{histone}/${histone}/g" ${dataDir}/${sampleID}_ChIP_seq_template_peak.json
			jbrowse add-track-json ${dataDir}/${sampleID}_ChIP_seq_template_peak.json --update

		elif [ $condition == "Control" -a $type == "bigwig" -a $experiment == "ChIP-seq" ] 
			then
			cp ${dataDir}/ChIP-seq/${sampleID}.bw ${jBrowseDir}/DATA/ChIP-seq
			sed "s/{sample}/${sampleID}/g"  ${dataDir}/ChIP_seq_control_template.json > ${dataDir}/${sampleID}_ChIP_seq_control_template.json 
			sed -i "s/{tissue}/${tissue}/g" ${dataDir}/${sampleID}_ChIP_seq_control_template.json 
			sed -i "s/{assemblyName}/${assemblyCode}/g" ${dataDir}/${sampleID}_ChIP_seq_control_template.json 
			sed -i "s/{histone}/${histone}/g" ${dataDir}/${sampleID}_ChIP_seq_control_template.json 
			jbrowse add-track-json ${dataDir}/${sampleID}_ChIP_seq_control_template.json  --update

		elif [ $condition == "Control" -a $type == "narrowPeak" -a $experiment == "ChIP-seq" ] 
			then
			cp ${dataDir}/ChIP-seq/${sampleID}.bigBed ${jBrowseDir}/DATA/ChIP-seq
			sed "s/{sample}/${sampleID}/g" ${dataDir}/ChIP_seq_control_peak_template.json > ${dataDir}/${sampleID}_ChIP_seq_control_peak_template.json 
			sed -i "s/{tissue}/${tissue}/g" ${dataDir}/${sampleID}_ChIP_seq_control_peak_template.json 
			sed -i "s/{assemblyName}/${assemblyCode}/g" ${dataDir}/${sampleID}_ChIP_seq_control_peak_template.json 
			sed -i "s/{histone}/${histone}/g" ${dataDir}/${sampleID}_ChIP_seq_control_peak_template.json 
			jbrowse add-track-json ${dataDir}/${sampleID}_ChIP_seq_control_peak_template.json  --update

		elif [ $condition != "-" -a $condition != "Control" -a $type == "bigwig" -a $experiment == "ChIP-seq" ] 
			then
			cp ${dataDir}/ChIP-seq/${sampleID}.bw ${jBrowseDir}/DATA/ChIP-seq
			sed "s/{sample}/${sampleID}/g" ${dataDir}/ChIP_seq_treatment_template.json > ${dataDir}/${sampleID}_ChIP_seq_treatment_template.json
			sed -i "s/{tissue}/${tissue}/g" ${dataDir}/${sampleID}_ChIP_seq_treatment_template.json
			sed -i "s/{assemblyName}/${assemblyCode}/g" ${dataDir}/${sampleID}_ChIP_seq_treatment_template.json 
			sed -i "s/{histone}/${histone}/g" ${dataDir}/${sampleID}_ChIP_seq_treatment_template.json 
			sed -i "s/{condition}/${condition}/g" ${dataDir}/${sampleID}_ChIP_seq_treatment_template.json 
			jbrowse add-track-json  ${dataDir}/${sampleID}_ChIP_seq_treatment_template.json   --update

		elif [ $condition != "-" -a $condition != "Control" -a $type == "narrowPeak" -a $experiment == "ChIP-seq" ] 
			then
			cp ${dataDir}/ChIP-seq/${sampleID}.bigBed ${jBrowseDir}/DATA/ChIP-seq
			sed "s/{sample}/${sampleID}/g" ${dataDir}/ChIP_seq_treatment_peak_template.json > ${dataDir}/${sampleID}_ChIP_seq_treatment_peak_template.json
			sed -i "s/{tissue}/${tissue}/g" ${dataDir}/${sampleID}_ChIP_seq_treatment_peak_template.json
			sed -i "s/{assemblyName}/${assemblyCode}/g" ${dataDir}/${sampleID}_ChIP_seq_treatment_peak_template.json
			sed -i "s/{histone}/${histone}/g" ${dataDir}/${sampleID}_ChIP_seq_treatment_peak_template.json
			sed -i "s/{condition}/${condition}/g" ${dataDir}/${sampleID}_ChIP_seq_treatment_peak_template.json
			jbrowse add-track-json  ${dataDir}/${sampleID}_ChIP_seq_treatment_peak_template.json   --update

		elif [ $condition == "Control" -a $type == "bigwig" -a $experiment == "ATAC-seq" ] 
			then
			cp ${dataDir}/ATAC-seq/${sampleID}.bw ${jBrowseDir}/DATA/ATAC-seq
			sed "s/{sample}/${sampleID}/g"  ${dataDir}/ATAC_seq_control_template.json > ${dataDir}/${sampleID}_ATAC_seq_control_template.json 
			sed -i "s/{tissue}/${tissue}/g" ${dataDir}/${sampleID}_ATAC_seq_control_template.json 
			sed -i "s/{assemblyName}/${assemblyCode}/g" ${dataDir}/${sampleID}_ATAC_seq_control_template.json 
			jbrowse add-track-json ${dataDir}/${sampleID}_ATAC_seq_control_template.json  --update

		elif [ $condition == "Control" -a $type == "narrowPeak" -a $experiment == "ATAC-seq" ] 
			then
			cp ${dataDir}/ATAC-seq/${sampleID}.bigBed ${jBrowseDir}/DATA/ATAC-seq
			sed "s/{sample}/${sampleID}/g" ${dataDir}/ATAC_seq_control_peak_template.json > ${dataDir}/${sampleID}_ATAC_seq_control_peak_template.json 
			sed -i "s/{tissue}/${tissue}/g" ${dataDir}/${sampleID}_ATAC_seq_control_peak_template.json 
			sed -i "s/{assemblyName}/${assemblyCode}/g" ${dataDir}/${sampleID}_ATAC_seq_control_peak_template.json 
			jbrowse add-track-json ${dataDir}/${sampleID}_ATAC_seq_control_peak_template.json  --update

		elif [ $condition != "-" -a $condition != "Control" -a $type == "bigwig" -a $experiment == "ATAC-seq" ] 
			then
			cp ${dataDir}/ATAC-seq/${sampleID}.bw ${jBrowseDir}/DATA/ATAC-seq
			sed "s/{sample}/${sampleID}/g" ${dataDir}/ATAC_seq_treatment_template.json > ${dataDir}/${sampleID}_ATAC_seq_treatment_template.json
			sed -i "s/{tissue}/${tissue}/g" ${dataDir}/${sampleID}_ATAC_seq_treatment_template.json
			sed -i "s/{assemblyName}/${assemblyCode}/g" ${dataDir}/${sampleID}_ATAC_seq_treatment_template.json 
			sed -i "s/{condition}/${condition}/g" ${dataDir}/${sampleID}_ATAC_seq_treatment_template.json 
			jbrowse add-track-json  ${dataDir}/${sampleID}_ATAC_seq_treatment_template.json   --update

		elif [ $condition != "-" -a $condition != "Control" -a $type == "narrowPeak" -a $experiment == "ATAC-seq" ] 
			then
			cp ${dataDir}/ATAC-seq/${sampleID}.bigBed ${jBrowseDir}/DATA/ATAC-seq
			sed "s/{sample}/${sampleID}/g" ${dataDir}/ATAC_seq_treatment_peak_template.json > ${dataDir}/${sampleID}_ATAC_seq_treatment_peak_template.json
			sed -i "s/{tissue}/${tissue}/g" ${dataDir}/${sampleID}_ATAC_seq_treatment_peak_template.json
			sed -i "s/{assemblyName}/${assemblyCode}/g" ${dataDir}/${sampleID}_ATAC_seq_treatment_peak_template.json
			sed -i "s/{condition}/${condition}/g" ${dataDir}/${sampleID}_ATAC_seq_treatment_peak_template.json
			jbrowse add-track-json  ${dataDir}/${sampleID}_ATAC_seq_treatment_peak_template.json   --update
		fi
done < ${dataDir}/data_ChIP-seq_ATAC-seq.tsv
