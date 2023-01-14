# DrosOmics
Scripts to create the DrosOmics genome browser: http://gonzalezlab.eu/drosomics

## Create environment
`conda create -n jbrowse2`

`conda activate jbrowse`
### Packages needed
`conda install gff3sort gffread tabix genometools samtools seqtk minimap2`

## Files
1. `createBrowser`: `jbrowse` instructions to load all data to DrosOmics
2. `data*tsv`: tab-delimited files with the name of the data uploaded into DrosOmics
3. `JSON-templates`: templates for DrosOmics to create the different type tracks
4. `ChIP-seq_scripts` and `RNA-seq_script`: scripts to process ChIP/ATAC/RNA-seq

## Citation

Coronado-Zamora M, Salces-Ortiz J, González J. 2022. DrosOmics: a comparative genomics browser to explore omics data in natural populations of *D. melanogaster*. BioRxiv:2022.07.22.501088. Available from: https://www.biorxiv.org/content/10.1101/2022.07.22.501088v2


