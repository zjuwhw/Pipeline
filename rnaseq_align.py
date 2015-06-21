#!/usr/local/bin/python
#author: zjuwhw
#2015-06-10
#for tophat, it seem there is a bug using both "tophat_report" and "--no-discordant". Some people said that tophat version2.0.9 and version 2.0.8b don't report the error.
USAGE='''rnaseq_align.py--- align pair-end strand-specific (or not) RNA-seq reads by Tophat and STAR
USAGE:
    python %s [--tool=#] [--strand-specific=#] [--tophat-genome-index=#] [--tophat-transcriptome-index=#] [--STAR-index=#] [--thread=#] [--tophat-other-options=#] [--STAR-other-options=#] [--chrom_size=#] [--dict=#] [--prefix=#] read1 read2

#input: quality and adapter trimmed fastq(.fastq, .fq, .fastq.gz, .fq.gz) files of read1 and read2
#output: $prefix_tophat folder, $prefix_STAR, $prefix_BamandSignal folder

#defaults:
--tool: "tophat", "STAR", "both". Default: "both"
--strand-specific: "stranded" or "unstranded". Default: "stranded". Now, this program can only deal with the dUTP method strand-specific pair-end RNA-seq, which is the method that Illumina stranded Trueseq takes.
--tophat-genome-index (bowtie2-index): "/d/database/hg38/DNAsequence_Ensembl/bowtie2index/Homo_sapiens.GRCh38.dna.primary_assembly".
--tophat-transcriptome-index: "/d/database/hg38/DNAsequence_Ensembl/TophatTranscriptomeIndex/Homo_sapiens.GRCh38.80.RNA".
--STAR-index: "/d/database/hg38/DNAsequence_Ensembl/STARindex/"
the tophat-genome-index, tophat-transcriptome-index, STAR-index can be built using another program in the Pipeline repo called "build-index.sh"
--thread: the number of threads is used. Default: 6.
--tophat-options: the options that tophat used.Default: -z0 -a 8 --min-intron-length 20 --max-intron-length 1000000 --read-edit-dist 4 --read-mismatches 4 -g 20 --no-discordant --no-mixed
The default options is what encode script used (https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/dnanexus/align-tophat-pe/src/align-tophat-pe.sh), but I do not include unmapped reads in the outputbam.
--STAR-options: the options that STAR used. Default:Default: --outSAMtype BAM Unsorted  --genomeLoad NoSharedMemory  --outFilterMultimapNmax 20  --alignSJoverhangMin 8 --alignSJDBoverhangMin 1  --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04  --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterType BySJout  --quantMode TranscriptomeSAM --sjdbScore 1 
The default options is what encode script used (https://github.com/ENCODE-DCC/long-rna-seq-pipeline/blob/master/dnanexus/align-star-pe/src/align-star-pe.sh), but I do not include unmapped reads in the outputbam.
--chrom_size: the chrom size of reference genome, which is used to convert bam to bigWig. Default: "/d/database/hg38/DNAsequence_Ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa.chrom.sizes"
--dict: the sequence dictionary of reference genome, which is used by picard. Default: "/d/database/hg38/DNAsequence_Ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa.dict"
the chrom_size file and dict file should be in the same directory with reference sequence and can be produced by another program in the Pipeline repo called "build-index.sh"
--prefix: the basename of read1 and trim the postfix "_1P.fq(fastq/fq.gz/fastq.gz)"

Note:
The tools, STAR, tophat, picard, samtools and UCSC utilities, are needed in $PATH
After alignment, the program will use the bam file to do:
1)SortSam by coodinate
2)add the RG(Read Group) to the bam header
3)Reorder the header of the bam to match with the sequence dictionary
4)SortSam by name
5)remove the middle files
6)convert bam to bigWig file to visualize in IGV.
7)creat soft link to bam and bigWig files in the $prefix_BamandSignal folder
'''

import os, getopt, sys, time, os.path, re

def tophat_stranded(read1, read2, TophatOpt, thread, prefix, TTI, TGI):
    print "### running tophat for stranded RNA-seq ..."
    cmd = "tophat --library-type fr-firststrand %s -p %d -o %s --transcriptome-index %s %s %s %s" % (TophatOpt, thread, prefix+"_tophat", TTI, TGI, read1, read2)
    os.system(cmd)
def tophat_unstranded(read1, read2, TophatOpt, thread, prefix, TTI, TGI):
    print "### running tophat for unstranded RNA-seq ..."
    cmd = "tophat --library-type fr-unstranded %s -p %d -o %s --transcriptome-index %s %s %s %s" % (TophatOpt, thread, prefix+"_tophat", TTI, TGI, read1, read2)
    os.system(cmd)
def STAR(read1, read2, STAROpt, thread, prefix, SI):
    print "### running STAR for stranded or unstranded RNA-seq ..."
    cmd1 = "mkdir -p "+prefix+"_STAR"
    os.system(cmd1)
    cmd2 = "STAR %s --runThreadN %d --genomeDir %s --readFilesIn %s %s --outFileNamePrefix %s" % (STAROpt, thread, SI, read1, read2, prefix+"_STAR/")
    os.system(cmd2)
def bam_to_bigWig_stranded(inputbam, prefix, chromsize):
    print "### convert stranded bam file to bigWig file ..."
    cmd1 = "STAR --runMode inputAlignmentsFromBAM --inputBAMfile %s --outWigType bedGraph --outWigStrand Stranded --outFileNamePrefix  %s/ --outWigNorm RPM" % (inputbam, os.path.dirname(inputbam))
    cmd2 = "bedGraphToBigWig %s/Signal.UniqueMultiple.str1.out.bg %s %s/%s_minusAll.bigWig" % (os.path.dirname(inputbam), chromsize, os.path.dirname(inputbam), prefix)
    cmd3 = "bedGraphToBigWig %s/Signal.Unique.str1.out.bg %s %s/%s_minusUniq.bigWig" % (os.path.dirname(inputbam), chromsize, os.path.dirname(inputbam), prefix)
    cmd4 = "bedGraphToBigWig %s/Signal.UniqueMultiple.str2.out.bg %s %s/%s_plusAll.bigWig" % (os.path.dirname(inputbam), chromsize, os.path.dirname(inputbam), prefix)
    cmd5 = "bedGraphToBigWig %s/Signal.Unique.str2.out.bg %s %s/%s_plusUniq.bigWig" % (os.path.dirname(inputbam), chromsize, os.path.dirname(inputbam), prefix)
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    os.system(cmd4)
    os.system(cmd5)
    os.remove("%s/Signal.UniqueMultiple.str1.out.bg" % os.path.dirname(inputbam))
    os.remove("%s/Signal.Unique.str1.out.bg" % os.path.dirname(inputbam))
    os.remove("%s/Signal.UniqueMultiple.str2.out.bg" % os.path.dirname(inputbam))
    os.remove("%s/Signal.Unique.str2.out.bg" % os.path.dirname(inputbam))
    
    if not os.path.exists(prefix+"_BamandSignal"):
        os.mkdir(prefix+"_BamandSignal")
    os.system("ln -s ../%s/%s_plusUniq.bigWig %s_BamandSignal/%s_plusUniq.bigWig" % (os.path.dirname(inputbam), prefix, prefix,os.path.dirname(inputbam).rstrip("/").split("/")[-1]))
    os.system("ln -s ../%s/%s_plusAll.bigWig %s_BamandSignal/%s_plusAll.bigWig" % (os.path.dirname(inputbam), prefix, prefix,os.path.dirname(inputbam).rstrip("/").split("/")[-1]))
    os.system("ln -s ../%s/%s_minusUniq.bigWig %s_BamandSignal/%s_minusUniq.bigWig" % (os.path.dirname(inputbam), prefix, prefix,os.path.dirname(inputbam).rstrip("/").split("/")[-1]))
    os.system("ln -s ../%s/%s_minusAll.bigWig %s_BamandSignal/%s_minusAll.bigWig" % (os.path.dirname(inputbam), prefix, prefix,os.path.dirname(inputbam).rstrip("/").split("/")[-1]))
def bam_to_bigWig_unstranded(inputbam, prefix, chromsize):
    print "### convert unstranded bam file to bigWig file ..."
    cmd1 = "STAR --runMode inputAlignmentsFromBAM --inputBAMfile %s --outWigType bedGraph --outWigStrand Unstranded --outFileNamePrefix  %s/ --outWigNorm RPM" % (inputbam, os.path.dirname(inputbam))
    cmd2 = "bedGraphToBigWig %s/Signal.UniqueMultiple.str1.out.bg %s %s/%s_all.bigWig" % (os.path.dirname(inputbam), chromsize, os.path.dirname(inputbam), prefix)
    cmd3 = "bedGraphToBigWig %s/Signal.Unique.str1.out.bg %s %s/%s_uniq.bigWig" % (os.path.dirname(inputbam), chromsize, os.path.dirname(inputbam), prefix)
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    os.remove("%s/Signal.UniqueMultiple.str1.out.bg" % os.path.dirname(inputbam))
    os.remove("%s/Signal.Unique.str1.out.bg" % os.path.dirname(inputbam))
    
    if not os.path.exists(prefix+"_BamandSignal"):
        os.mkdir(prefix+"_BamandSignal")
    os.system("ln -s ../%s/%s_uniq.bigWig %s_BamandSignal/%s_uniq.bigWig" % (os.path.dirname(inputbam), prefix, prefix, os.path.dirname(inputbam).rstrip("/").split("/")[-1]))
    os.system("ln -s ../%s/%s_all.bigWig %s_BamandSignal/%s_all.bigWig" % (os.path.dirname(inputbam), prefix, prefix, os.path.dirname(inputbam).rstrip("/").split("/")[-1]))
def indexbam(inputbam):
    cmd = "samtools index %s" % inputbam
    os.system(cmd)
def after_bam(inputbam, dictfile, strandness, prefix, chromsize):
    #the picard is usually used by the long command "java -jar [path to picard.jar]"
    #In this program, I write a script to use "picard" instead of long command:
    ##!/bin/bash
    #java -Xmx2g -jar /c/wanghw/software/rnaseq/picard-tools-1.133/picard.jar $@
    
    coordsortedbam = inputbam.replace(".bam", ".CoordSorted.tmp.bam")
    cmd1 = "samtools sort %s %s" % (inputbam, coordsortedbam.replace(".bam",""))
    os.system(cmd1)
    indexbam(coordsortedbam)
    os.remove(inputbam)
    
    RGroupbam = inputbam.replace(".bam", ".CoordSorted.RG.bam")
    cmd2 = "picard AddOrReplaceReadGroups I=%s O=%s ID=Wanglab-%s LB=%s PL=ILLUMINA PU=a SM=%s" % (coordsortedbam, RGroupbam, prefix, strandness, prefix)
    os.system(cmd2)
    indexbam(RGroupbam)
    os.remove(coordsortedbam)
    os.remove(coordsortedbam+".bai")
    
    RSamheader = inputbam.replace(".bam", ".CoordSorted.bam")
    cmd3 = "picard ReorderSam I=%s O=%s R=%s" % (RGroupbam, RSamheader, dictfile.replace(".dict",""))
    os.system(cmd3)
    indexbam(RSamheader)
    os.remove(RGroupbam)
    os.remove(RGroupbam+".bai")
    
    namesortedbam = inputbam.replace(".bam", ".NameSorted.bam")
    cmd4 = "samtools sort -n %s %s" % (RSamheader, namesortedbam.replace(".bam",""))
    os.system(cmd4)   
    
    if strandness == "stranded":
        bam_to_bigWig_stranded(RSamheader, prefix, chromsize)
    elif strandness == "unstranded":
        bam_to_bigWig_unstranded(RSamheader, prefix, chromsize)    
    
    if not os.path.exists(prefix+"_BamandSignal"):
        os.mkdir(prefix+"_BamandSignal")
    os.system("ln -s ../%s/%s %s_BamandSignal/%s.CoordSorted.bam" % (os.path.dirname(inputbam), os.path.basename(inputbam).replace(".bam",".CoordSorted.bam"),prefix,os.path.dirname(inputbam).rstrip("/").split("/")[-1]))
    os.system("ln -s ../%s/%s %s_BamandSignal/%s.CoordSorted.bam.bai" % (os.path.dirname(inputbam), os.path.basename(inputbam).replace(".bam",".CoordSorted.bam.bai"),prefix,os.path.dirname(inputbam).rstrip("/").split("/")[-1]))
    os.system("ln -s ../%s/%s %s_BamandSignal/%s.NameSorted.bam" % (os.path.dirname(inputbam), os.path.basename(inputbam).replace(".bam",".NameSorted.bam"),prefix,os.path.dirname(inputbam).rstrip("/").split("/")[-1]))
    
if __name__ == '__main__':
    
    starttime=time.time()
    if len(sys.argv) < 2:
        print USAGE % sys.argv[0]
        sys.exit(1)
    
    opts, args = getopt.getopt(sys.argv[1:], "", ["tool=","strand-specific=","tophat-genome-index=","tophat-transcriptome-index=","STAR-index=","thread=","tophat-other-options=", "STAR-other-options=", "chrom-size=","dict=", "prefix="])
    
    if len(args) != 2:
        print "This program is only for pair-end RNA-seq, not single-end RNA-seq"
        sys.exit(1)
    
    read1 = args[0]
    read2 = args[1]
    
    # defaults
    tools = ["both", "tophat", "STAR"]
    tool = "both"
    strandnesses = ["stranded", "unstranded"]
    strandness = "stranded"
    TGI = "/d/database/hg38/DNAsequence_Ensembl/bowtie2index/Homo_sapiens.GRCh38.dna.primary_assembly"
    TTI = "/d/database/hg38/DNAsequence_Ensembl/TophatTranscriptomeIndex/Homo_sapiens.GRCh38.80.RNA"
    SI = "/d/database/hg38/DNAsequence_Ensembl/STARindex/"
    thread = 6
    TophatOpt = "-z0 -a 8 --min-intron-length 20 --max-intron-length 1000000 --read-edit-dist 4 --read-mismatches 4 -g 20 --no-discordant --no-mixed "
    STAROpt = "--outSAMtype BAM Unsorted  --genomeLoad NoSharedMemory  --outFilterMultimapNmax 20  --alignSJoverhangMin 8 --alignSJDBoverhangMin 1  --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04  --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterType BySJout  --quantMode TranscriptomeSAM --sjdbScore 1 "
    chromsize = "/d/database/hg38/DNAsequence_Ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa.chrom.sizes"
    dictfile = "/d/database/hg38/DNAsequence_Ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa.dict"
    regular = "_[12]\.(fastq|fq)(\.gz)?$"
    if (not re.search(regular, os.path.basename(read1))) or (not re.search(regular, os.path.basename(read2))):
        print '''the input fastq file name must obey the regular expression "_[12]\.(fastq|fq)(\.gz)?$"'''
        sys.exit(1)
    prefix = re.sub(regular, "", os.path.basename(read1))
    if os.path.basename(read1).endswith(".gz"):
        STAROpt += " --readFilesCommand zcat "
    
    for o, a in opts:
        if o == "--tool":
            tool = a
            if tool not in tools:
                print "Only three tool options are avalible, tophat, STAR, both "
                sys.exit(1)
        if o == "--strand-specific":
            strandness = a
            if strandness not in strandnesses:
                print "Only two strand-specific options are avalible, stranded or unstranded"
                sys.exit(1)
        if o == "--tophat-genome-index":
            TGI = a
        if o == "--tophat-transcriptome-index":
            TTI = a
        if o == "--STAR-index" :
            SI = a
        if o == "--thread" :
            thread = int(a)
        if o == "--tophat-other-options":
            TophatOpt = a
        if o == "--STAR-other-options" :
            STAROpt = a
        if o == "--chrom-size":
            chromsize = a
        if o == "--dict":
            dictfile = a
        if o == "--prefix" :
            prefix = a
    
    if tool == "both":
        if strandness == "stranded":
            tophat_stranded(read1, read2, TophatOpt, thread, prefix, TTI, TGI)
            STAR(read1, read2, STAROpt, thread, prefix, SI)
        elif strandness == "unstranded":
            tophat_unstranded(read1, read2, TophatOpt, thread, prefix, TTI, TGI)
            STAR(read1, read2, STAROpt, thread, prefix, SI)
        after_bam(prefix+"_tophat/accepted_hits.bam", dictfile, strandness, prefix, chromsize)
        after_bam(prefix+"_STAR/Aligned.out.bam", dictfile, strandness, prefix, chromsize)
    elif tool == "tophat":
        if strandness == "stranded":
            tophat_stranded(read1, read2, TophatOpt, thread, prefix, TTI, TGI)
        elif strandness == "unstranded":
            tophat_unstranded(read1, read2, TophatOpt, thread, prefix, TTI, TGI)
        after_bam(prefix+"_tophat/accepted_hits.bam", dictfile, strandness, prefix, chromsize)
    elif tool == "STAR":
        if strandness == "stranded":
            STAR(read1, read2, STAROpt, thread, prefix, SI)
        elif strandness == "unstranded":
            STAR(read1, read2, STAROpt, thread, prefix, SI)
        after_bam(prefix+"_STAR/Aligned.out.bam", dictfile, strandness, prefix, chromsize)
    
    endtime = time.time()
    print "it takes %.3f minutes or %.3f seconds" % ((endtime-starttime)/60, endtime -starttime)