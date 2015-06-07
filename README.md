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
*chipseq_pipeline.py for alignment using the sra or fastq files, and output bam file with only the unique mapped and non duplicated reads

*chipseq_peakcalling.py for peak calling using unique mapped and non duplicated reads. Two tools are avaliable, MACS14 and Hpeak.


*chipseq_peak_overlap_venn.py for cacluate the overlap number of different ChIP-seq peaks and then do the Venn plot

*chipseq_peak2gene.py for find the target gene of a ChIPed enrichment region (peak). Three types to associate the peak with gene are avaliable, that is "peak2gene", "gene2peak", "peakAroundgene"

*chipseq_extract_signal_from_bigwig.py for calculate the signal value for a genome region and then to do aggregation plot or heatmap plot

*chipseq_motif_discovery.py for motif discovery

### Microarray analysis:
####required tools:

