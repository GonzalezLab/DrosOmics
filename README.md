# DrosOmics
Docker image and scripts to run and to create the DrosOmics genome browser: http://gonzalezlab.eu/drosomics

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

## Docker image
The docker image containing DrosOmics is in the `drosomics.tar` file.

## Citation

Coronado-Zamora M, Salces-Ortiz J, González J. 2023. DrosOmics: a browser to explore -omics variation across high-quality reference genomes from natural populations of Drosophila melanogaster. Mol Biol Evol 40: msad075. 10.1093/molbev/msad075
