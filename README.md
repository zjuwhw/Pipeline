Pipeline
========

This is a repository to store the pipeline to analysis affy microarray, single-end ChIP-seq and pair-end RNA-seq in Wang lab.
### ChIP-seq analysis:
####required tools:
*alignment: bwa, bowtie2

*peak calling: MACS, hpeak

*motif discovery: meme-chip, homer, amd

*plot: VennDiagram, ggplot2 in R software

*other tools: bedtools, samtools, UCSC Jim Kent utility

####pipeline:
*chipseq_pipeline.py for alignment using the sra or fastq files, and output bam file with only the uniq mapped and no duplicated reads

*extract_signal_from_bigwig.py for calculate the signal value for a genome region and then to do aggregation plot or heatmap plot

*motif_discovery.py for motif discovery

*peak_overlap_venn.py for cacluate the overlap number of different ChIP-seq peaks and then do the Venn plot

### Microarray analysis:
####required tools:

