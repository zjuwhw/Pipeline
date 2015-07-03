#!/usr/local/bin/python
#author: zjuwhw
#email: zju.whw@gmail.com
#2015-07-02
USAGE='''rnaseq_quant.py--- quantification of RNA-seq bam files using RSEM (for transcritome bam) or featureCount (for genome bam)
USAGE:
    python %s [--tool=#] [--rsemindex=#] [--gtf=#] [--thread=#] [--strand-specific=#] [--prefix=#] inputbam

#input: a aligned bam files. For RSEM, the inputbam is the one aligned to the transcriptome produced by STAR. For featureCount, the inputbam is a name sorted one aligend genome and transcriptome by STAR or tophat.
#output: $prefix_rsem folder and/or $prefix_count folder

#defaults:
--tool: "rsem" or "featurecount". Default: "rsem"
--rsemindex: the path to the rsem index prefix, it can be built using rsem-prepare-reference. Default:/d/database/hg38/DNAsequence_Ensembl/rsemindex/resm.
--gtf: the gene annotation. Default: /d/database/hg38/GeneAnn_Ensembl/Homo_sapiens.GRCh38.80.gtf
--thread: the number of threads is used. Default: 8.
--strand-specific: "stranded" or "unstranded". Default: "stranded". Now, this program can only deal with the dUTP method strand-specific pair-end RNA-seq, which is the method that Illumina stranded Trueseq takes.
--prefix: the basename of inputbam

Note:
RSEM and featureCount should be built in the $PATH
'''

import os, getopt, sys, time, os.path, re

def mkdir(dirctory):
    if not os.path.exists(dirctory):
        os.mkdir(dirctory)

def rsem(inputbam, RsemOpt, rsemindex, prefix):
    starttime_rsem=time.time()
    print "#########rsem is beginning ..."
    mkdir(prefix+"_rsem")
    cmd = "rsem-calculate-expression --bam %s %s %s %s" % (RsemOpt, inputbam, rsemindex, prefix+"_rsem/"+prefix)
    os.system(cmd)
    endtime_rsem = time.time()
    print "the rsem step takes %.3f minutes or %.3f seconds" % ((endtime_rsem - starttime_rsem)/60, endtime_rsem -starttime_rsem)
def featurecount(inputbam, feaOpt, gtf, prefix):
    starttime_featurecount=time.time()
    print "#########featurecount is beginning ..."
    mkdir(prefix+"_count")
    cmd1 = "featureCounts -p %s -a %s -o %s_count/%s_gene.fragmentcount %s" % (feaOpt, gtf, prefix, prefix, inputbam)
    cmd2 = "featureCounts %s -a %s -o %s_count/%s_gene.readcount %s" % (feaOpt, gtf, prefix, prefix, inputbam)
    os.system(cmd1)
    os.system(cmd2)
    endtime_featurecount = time.time()
    print "the featurecount step takes %.3f minutes or %.3f seconds" % ((endtime_featurecount - starttime_featurecount)/60, endtime_featurecount -starttime_featurecount)
    

if __name__ == '__main__':
    
    starttime=time.time()
    if len(sys.argv) < 2:
        print USAGE % sys.argv[0]
        sys.exit(1)
        
    
    opts, args = getopt.getopt(sys.argv[1:], "", ["tool=", "gtf=", "rsemindex=","strand-specific=", "thread=", "prefix="])
    
    inputbam = args[0]
    
    #default:
    tools = ["rsem", "featurecount"]
    tool = "rsem"
    strandnesses = ["stranded", "unstranded"]
    strandness = "stranded"
    rsemindex = "/d/database/hg38/DNAsequence_Ensembl/rsemindex/resm"
    gtf = "/d/database/hg38/GeneAnn_Ensembl/Homo_sapiens.GRCh38.80.gtf"
    thread = 8
    regular = "\.bam$"
    if not re.search(regular, os.path.basename(inputbam)):
        print '''the input is a bam file end with ".bam"'''
        sys.exit(1)
    prefix = re.sub(regular, "", os.path.basename(inputbam))
    
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
        if o == "--thread" :
            thread = int(a)
        if o == "--gtf" :
            gtf = a
        if o == "--rsemindex" :
            rsemindex = a
        if o == "--prefix":
            prefix = a
    
    
    
    if tool == "rsem":
        if strandness == "stranded":
            RsemOpt = " -p %d --no-bam-output --paired-end --forward-prob 0 "  % thread
        elif strandness == "unstranded":
            RsemOpt = " -p %d --no-bam-output --paired-end --forward-prob 0.5 "  % thread
        rsem(inputbam, RsemOpt, rsemindex, prefix)
    elif tool == "featurecount":
        if strandness == "stranded":
            feaOpt = " -T %d -s 2 " % thread
        elif strandness == "unstranded":
            feaOpt = " -T %d -s 0 " % thread
        featurecount(inputbam, feaOpt, gtf, prefix)
    
    endtime = time.time()
    print "it takes %.3f minutes or %.3f seconds" % ((endtime-starttime)/60, endtime -starttime)