Pipeline
========

This is a repository to store the pipeline to analysis affy microarray, single-end ChIP-seq and pair-end RNA-seq in Wang lab.
### single-end ChIP-seq analysis:
####required tools:
*alignment: [bwa](http://bio-bwa.sourceforge.net/), [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

*peak calling: [MACS](http://liulab.dfci.harvard.edu/MACS/), [hpeak](http://www.sph.umich.edu/csg/qin/HPeak/)

*motif discovery: [meme-chip](http://meme.nbcr.net/meme/intro.html), [homer](http://homer.salk.edu/homer/), [amd](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0024576)

*plot: [VennDiagram](http://cran.r-project.org/web/packages/VennDiagram/index.html), [ggplot2](http://ggplot2.org/) in R software

*other tools: [bedtools](http://bedtools.readthedocs.org/en/latest/), [samtools](http://samtools.sourceforge.net/), [UCSC Jim Kent utility](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/)

####pipeline:
*chipseq_pipeline.py --- alignment using the sra or fastq files, and output bam file with only the unique mapped and non duplicated reads

*chipseq_peakcalling.py --- peak calling using unique mapped and non duplicated reads. Two tools are avaliable, MACS14 and Hpeak.

*chipseq_peak_overlap_venn.py --- cacluate the overlap number of different ChIP-seq peaks and then do the Venn plot

*chipseq_peak2gene.py --- find the target gene of a ChIPed enrichment region (peak). Three types to associate the peak with gene are avaliable, that is "peak2gene", "gene2peak", "peakAroundgene"

*chipseq_extract_signal_from_bigwig.py for calculate the signal value for a genome region and then to do aggregation plot or heatmap plot

*chipseq_motif_discovery.py --- motif discovery using meme-chip, homer or amd

### Microarray analysis:
####required tools:
*R bioconductor package [affy](http://www.bioconductor.org/packages/release/bioc/html/affy.html)
for "hgu133a","hgu133a2","hgu133b","hgu133plus2","hgu219","hgu95a","hgu95av2","hgu95b","hgu95c","hgu95d","hgu95e","u133aaofav2"

*R bioconductor package [oligo](http://www.bioconductor.org/packages/release/bioc/html/oligo.html)
for "huex10st","hugene10st","hugene11st","hugene20st","hugene21st"

####pipeline:
*affy_build_annotation.py --- download annotation files for affy microarray probe id, using the R bioconductor annotation db package 

*affy_ExonOrGene_build_annotation.py --- download annotation files for affy Exon Or Gene microarray transcriptcluster id, using the R bioconductor annotation db package

*affy_array_pipeline.py --- affy microarray pipeline using R to do rma, mas5.0 and/or not customCDF normalization

*affy_array_pipeline.py --- affy microarray pipeline using R to do rma, mas5.0 and/or not customCDF normalization

###pair-end RNA-seq analysis:
####required tools:
*fastq QC: [fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

*trimming: [trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)

*alignment: [STAR](https://github.com/alexdobin/STAR) [tophat](http://tophat.cbcb.umd.edu/)

*bam QC:
