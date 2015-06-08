USAGE='''rnaseq_QCandTrimming.py--- Qualtiy Control and Trimming for Illumina PAIRED-end RNAseq fastq files.
Quality Control for fastq before and after trimming using FastQC,
Qualtiy and adapter trimming using trimmomatic.

USAGE:
    python %s [--trimmomatic_option=#][--adapter=#] [--outDir=#] read1_fastq read2_fastq

#input: fastq(.fastq, .fq, .fastq.gz, .fq.gz) files of read1 and read2
#output: 4 fastqc.html files and fastq files ofpaired(P)/unpaired(U) read1 and read2

#defaults:
#--trimmomatic_option: ILLUMINACLIP:$adapter:2:30:10 SLIDINGWINDOW:4:15 TRAILING:10 MINLEN:20, which means
1)Remove Illumina adapters provided in the adapter fasta file. Initially Trimmomatic will look for seed matches (16 bases) allowing maximally 2 mismatches. These seeds will be extended and clipped if in the case of paired end reads a score of 30 is reached (about 50 bases), or in the case of single end reads a score of 10, (about 17 bases)
2)Remove trailing low quality or N bases (below quality 10)
3)scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15.
4)Drop reads which are less than 20 bases long after these steps.
#--adapter: the path of adapter fasta file, trimmomatic provides "TruSeq2-PE.fa" for GAII machines and "TruSeq3-PE.fa" or "TruSeq3-PE-2.fa" for HiSeq and MiSeq machines. Default is "/c/wanghw/software/rnaseq/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa".
#outDir: current directory

Note:
The tools, fastqc and trimmomatic, are needed in $PATH
'''

import os, getopt, sys, time, os.path

def fastqc(fastqfile, outdir):
    cmd = "fastqc -o %s %s" % (outdir, fastqfile)
    os.system(cmd)

def trim(read1, read2, outdir, option):
    read1P_name=os.path.basename(read1).replace("1.fastq","1P.fastq").replace("1.fastq.gz","1P.fastq.gz").replace("1.fq.gz", "1P.fq.gz").replace("1.fq", "1P.fq")
    read1U_name=os.path.basename(read1).replace("1.fastq","1U.fastq").replace("1.fastq.gz","1U.fastq.gz").replace("1.fq.gz", "1U.fq.gz").replace("1.fq", "1U.fq")
    read2P_name=os.path.basename(read2).replace("2.fastq","2P.fastq").replace("2.fastq.gz","2P.fastq.gz").replace("2.fq.gz", "2P.fq.gz").replace("2.fq", "2P.fq")
    read2U_name=os.path.basename(read2).replace("2.fastq","2U.fastq").replace("2.fastq.gz","2U.fastq.gz").replace("2.fq.gz", "2U.fq.gz").replace("2.fq", "2U.fq")
    
    cmd = "trimmomatic PE %s %s %s %s %s %s %s" % (read1, read2, outdir.rstrip("/")+"/"+read1P_name, outdir.rstrip("/")+"/"+read1U_name, outdir.rstrip("/")+"/"+read2P_name, outdir.rstrip("/")+"/"+read2U_name, option)
    print "########Trimmomatic command:" + cmd
    os.system(cmd)
    return [read1P_name, read2P_name]
    
if __name__ == '__main__':
    
    starttime=time.time()
    if len(sys.argv) < 2:
        print USAGE % sys.argv[0]
        sys.exit(1)
        
    opts, args = getopt.getopt(sys.argv[1:], "", ["trimmomatic_option=","adapter=","outDir="])
    
    if len(args) != 2:
        print "This program is only for pair-end RNA-seq, not single-end RNA-seq"
        sys.exit(1)
    
    read1 = args[0]
    read2 = args[1]
    
    # defaults
    outDir = "./"
    adapter="/c/wanghw/software/rnaseq/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa"
    option="ILLUMINACLIP:%s:2:30:10 SLIDINGWINDOW:4:15 TRAILING:10 MINLEN:20" % adapter
    for o,a in opts:
        if o == "--adapter":
            adapter = a
            option = "ILLUMINACLIP:%s:2:30:10 SLIDINGWINDOW:4:15 TRAILING:10 MINLEN:20" % a
        elif o == "--trimmomatic_option":
            option = a
        elif o == "outDir":
            outDir = a
    
    
    fastqc(read1, outDir)
    fastqc(read2, outDir)
    read1P, read2P = trim(read1, read2, outDir, option)
    print read1P, read2P
    fastqc(read1P, outDir)
    fastqc(read2P, outDir)

    endtime = time.time()
    print "it takes %.3f secondes or %.3f minutes" % (endtime -starttime, (endtime-starttime)/60)