#!/usr/local/bin/python
#author: zjuwhw
#2015-06-2

USAGE='''rnaseq_bamQC.py--- quality control for aligned RNA-seq bam files
USAGE:
    python %s [--tool=#] [--strand-specific=#] [--prefix=#] [--generefFlat=#] [--genegtf=#] [--genebed=#] [--rRNA_interval_list=#] [--dnaref=#] inputbam

#input: input aligned RNA-seq bam file
#output: output four folders, including ${prefix}_samstat, ${prefix}_rnaseqmetrics, ${prefix}_rseqc, ${prefix}_rnaseqc

#defaults:
--tool: "samstat", "picard-rnaseqmetrics", "rseqc", "rnaseqc", "all"
--strand-specific: "stranded" or "unstranded". Default: "stranded". Now, this program can only deal with the dUTP method strand-specific pair-end RNA-seq, which is the method that Illumina stranded Trueseq takes.
--generefFlat: "/d/database/hg38/GeneAnn_Ensembl/Homo_sapiens.GRCh38.80.refFlat"
--genegtf: "/d/database/hg38/GeneAnn_Ensembl/Homo_sapiens.GRCh38.80.modified_RNASeQC.gtf"
--genebed: "/d/database/hg38/GeneAnn_Ensembl/Homo_sapiens.GRCh38.80.bed"
--rRNA_interval_list: "/d/database/hg38/GeneAnn_Ensembl/Homo_sapiens.GRCh38.80.rRNA.interval.list"
--dnaref: "/d/database/hg38/DNAsequence_Ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
--prefix: basename of inputbam file

Note:
The tools, SAMSTAT, picard, RSeQC and RNA-SeQC, are needed in $PATH
'''

import os, getopt, sys, time, os.path, re
def mkdir(foldername):
    if not os.path.exists(foldername):
        os.mkdir(foldername)
def samstat(inputbam, prefix):
    print "#########samstat is beginning ..."
    mkdir(prefix+"_samstat")
    os.system("samstat %s" % inputbam)
    os.system("mv %s %s" % (inputbam+".samstat.html", prefix+"_samstat"))
    
def picard(inputbam, prefix, generefFlat, rRNA_interval_list, picard_strandness):
    print "#########picard-collectrnaseqmetrics is beginning ..."
    mkdir(prefix+"_rnaseqmetrics")
    cmd = "picard CollectRnaSeqMetrics REF_FLAT=%s RIBOSOMAL_INTERVALS=%s STRAND_SPECIFICITY=%s MINIMUM_LENGTH=100 CHART_OUTPUT=%s INPUT=%s OUTPUT=%s" % (generefFlat, rRNA_interval_list, picard_strandness, prefix+".picard_output.pdf", inputbam, prefix+".picard_output.txt")
    os.system(cmd)
    os.system("mv %s %s" % (prefix+".picard_output.pdf", prefix+"_rnaseqmetrics"))
    os.system("mv %s %s" % (prefix+".picard_output.txt", prefix+"_rnaseqmetrics"))

def rnaseqc(inputbam, prefix, dnaref, genegtf, rRNA_interval_list):
    print "#########rnaseqqc is beginning ..."
    cmd = '''rnaseqc -r %s -s "XXX|%s|XXX" -t %s -o %s -rRNA %s''' % (dnaref, inputbam, genegtf, prefix+"_rnaseqc", rRNA_interval_list)
    os.system(cmd)

def rseqc(inputbam, prefix, genebed, rseqc_strandness):
    print "#########rseqc is beginning ..."
    mkdir(prefix+"_rseqc")
    
    print "#########rseqc-inner_distance.py is beginning ..."
    cmd1 = "geneBody_coverage.py -r %s -i %s -o %s" % (genebed, inputbam, prefix+"_rseqc/"+prefix)
    os.system(cmd1)
    
    print "#########rseqc-geneBody_coverage.py is beginning ..."
    cmd2 = "inner_distance.py -r %s -i %s -o %s" % (genebed, inputbam, prefix+"_rseqc/"+prefix)
    os.system(cmd2)
    
    print "#########rseqc-junction_annotation.py is beginning ..."
    cmd3 = "junction_annotation.py -r %s -i %s -o %s" % (genebed, inputbam, prefix+"_rseqc/"+prefix)
    os.system(cmd3)
    
    print "#########rseqc-junction_saturation.py is beginning ..."
    cmd4 = "junction_saturation.py -r %s -i %s -o %s" % (genebed, inputbam, prefix+"_rseqc/"+prefix)
    os.system(cmd4)
    
    print "#########rseqc-read_distribution.py is beginning ..."
    cmd5 = "read_distribution.py -r %s -i %s > %s" % (genebed, inputbam, prefix+"_rseqc/"+prefix + ".read_distribution.txt")
    os.system(cmd5)
    
    print "#########rseqc-read_duplication.py is beginning ..."
    cmd6 = "read_duplication.py -i %s -o %s" % (inputbam, prefix+"_rseqc/"+prefix + ".read_duplication")
    os.system(cmd6)
    
    print "#########rseqc-read_GC.py is beginning ..."
    cmd7 = "read_GC.py -i %s -o %s" % (inputbam, prefix+"_rseqc/"+prefix + ".read_GC")
    os.system(cmd7)
    
    print "#########rseqc-RPKM_saturation.py is beginning ..."
    cmd8 = 'RPKM_saturation.py -d "%s" -r %s -i %s -o %s' % (rseqc_strandness, genebed, inputbam, prefix+"_rseqc/"+prefix + ".RPKM_saturation")
    os.system(cmd8)

    
if __name__ == '__main__':
    starttime=time.time()
    if len(sys.argv) < 2:
        print USAGE % sys.argv[0]
        sys.exit(1)
    
    opts, args = getopt.getopt(sys.argv[1:], "", ["tool=", "strand-specific=", "prefix=", "generefFlat=", "genegtf=", "rRNA_interval_list=", "dnaref=", "genebed="])
    inputbam = args[0]
    
    #default:
    tools = ["samstat", "picard-rnaseqmetrics", "rseqc", "rnaseqc", "all"]
    tool = "all"
    strandnesses = ["stranded", "unstranded"]
    strandness = "stranded"
    prefix = os.path.basename(inputbam).replace(".bam", "")
    generefFlat = "/d/database/hg38/GeneAnn_Ensembl/Homo_sapiens.GRCh38.80.refFlat"
    genegtf = "/d/database/hg38/GeneAnn_Ensembl/Homo_sapiens.GRCh38.80.modified_RNASeQC.gtf"
    genebed= "/d/database/hg38/GeneAnn_Ensembl/Homo_sapiens.GRCh38.80.bed"
    rRNA_interval_list = "/d/database/hg38/GeneAnn_Ensembl/Homo_sapiens.GRCh38.80.rRNA.interval.list"
    dnaref = "/d/database/hg38/DNAsequence_Ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    
    
    
    
    for o, a in opts:
        if o == "--tool":
            tool = a
            if tool not in tools:
                print 'Only five tool options are avalible, "samstat", "picard-rnaseqmetrics", "rseqc", "rnaseqc", "all"'
                sys.exit(1)
        if o == "--strand-specific":
            strandness = a
            if strandness not in strandnesses:
                print "Only two strand-specific options are avalible, stranded or unstranded"
                sys.exit(1)
        if o == "--prefix":
            prefix = a
        if o == "--generefFlat":
            generefFlat = a
        if o =="--genegtf":
            genegtf = a
        if o == "--rRNA_interval_list":
            rRNA_interval_list = a
        if o == "--dnaref":
            dnaref = a
        if o == "--genebed":
            genebed = a
    
    picard_strandnesses = {"unstranded": "NONE", "stranded": "SECOND_READ_TRANSCRIPTION_STRAND"}
    picard_strandness = picard_strandnesses[strandness]
    rseqc_strandnesses = {"unstranded": "none", "stranded": "1+-,1-+,2++,2--"}
    rseqc_strandness = rseqc_strandnesses[strandness]
    
    switch_tools = {"samstat": False, "picard-rnaseqmetrics": False, "rseqc": False, "rnaseqc": False}
    
    if tool == "all":
        switch_tools = {"samstat": True, "picard-rnaseqmetrics": True, "rseqc": True, "rnaseqc": True}
    else:
        switch_tools[tool] = True
    
    if switch_tools["samstat"]:
        samstat(inputbam, prefix)
    if switch_tools["picard-rnaseqmetrics"]:
        picard(inputbam, prefix, generefFlat, rRNA_interval_list, picard_strandness)
    if switch_tools["rnaseqc"]:
        rnaseqc(inputbam, prefix, dnaref, genegtf, rRNA_interval_list)
    if switch_tools["rseqc"]:
        rseqc(inputbam, prefix, genebed, rseqc_strandness)  
        