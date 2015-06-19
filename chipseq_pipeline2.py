##author: zjuwhw
##time: 2015-06-19
USAGE='''chipseq_pipeline2.py --- This program is used for single-end ChIP-seq data preprocess and alignment by bwa software;
usage:
    python %s [--refdna=#] [--prefix=#] inputdata(fastq or sra)

#input: fastq/fq(.gz) or sra file, it must obey the regular expression "\.((fastq|fq)(\.gz)?|sra)$"
#output: ${prefix}.bam, ${prefix}.bam.bai, ${prefix}.log files

#defaults:
--refdna: DNA reference sequence fasta file. Default: "/d/database/hg19/bwaindex/hg19.fa".
--prefix: Default: the basename of inputdata.

Note:
The program "bwa","samtools" must be installed and avaliable in the $PATH.
BWA index of DNA reference sequence fasta must be built before.
'''

import os, getopt, sys, time, os.path, re

def ossystemresult(command):
    fp=os.popen(command,"r")
    return fp.read().rstrip()
def sra2fastq(srafile, prefix):
    cmd = "fastq-dump -Z %s > %s" % (srafile, prefix+".fastq")
    os.system(cmd)
    return prefix+".fastq"    
def fastq2bam(fastqfile, prefix, refdna):
    cmd1="bwa aln -t 4 %s %s > %s" % (refdna, fastqfile, prefix+".sai")
    cmd2="bwa samse %s %s %s > %s" % (refdna, prefix+".sai", fastqfile, prefix+".sam")
    cmd3="samtools view -bS %s > %s" % (prefix+".sam", prefix+"_unsorted.bam")
    cmd4="samtools sort %s %s" % (prefix+"_unsorted.bam", prefix)
    cmd5="rm -f *sam *sai *_unsorted.bam"
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    os.system(cmd4)
    os.system(cmd5)
    return prefix + ".bam"
def bam2uniqbam(bamfile):
    cmd = "samtools view -bq 1 %s > %s" % (bamfile, bamfile.replace(".bam","_uniq.bam"))
    os.system(cmd)
    return bamfile.replace(".bam","_uniq.bam")
def bam2nodupbam(bamfile):
    cmd="samtools rmdup -s %s %s" % (bamfile, bamfile.replace(".bam", "_nodup.bam"))
    os.system(cmd)
    return bamfile.replace(".bam", "_nodup.bam")
def bamindex(bamfile):
    cmd = "samtools index %s" % bamfile
    os.system(cmd)
def bamcount(bamfile):
    cmd = "samtools view %s | wc -l " % bamfile
    return ossystemresult(cmd)
        
if __name__ == '__main__':
    
    starttime=time.time()
    if len(sys.argv) < 2:
        print USAGE % sys.argv[0]
        sys.exit(1)
    
    opts, args = getopt.getopt(sys.argv[1:], "", ["refdna=", "prefix="])
    inputfile = args[0]
    
    regular = r"\.((fastq|fq)(\.gz)?|sra)$"
    if not re.search(regular, inputfile):
        print '''the input file must obey the regular expression "\.((fastq|fq)(\.gz)?|sra)$" '''
        sys.exit(1)
    postfix = re.search(regular, inputfile).group(0)
    
    #default:
    refdna = "/d/database/hg19/bwaindex/hg19.fa"
    prefix = re.sub(regular, "", os.path.basename(inputfile))
    
    for o, a in opts:
        if o == "--refdna":
            refdna = a
        if o == "--prefix":
            prefix = a
    
    logfile = open(prefix+".log", "w")
    print >> logfile, "inputfile name:\t%s" % (inputfile)
    print >> logfile, "inputfile type:\t%s" % (postfix)
    
    if postfix == ".sra":
        inputfile = sra2fastq(inputfile, prefix)
    bamfile = fastq2bam(inputfile, prefix, refdna)
    nodupbamfile = bam2nodupbam(bamfile)
    nodupuniqbamfile=bam2uniqbam(nodupbamfile)    
    bamindex(nodupuniqbamfile)
    
    print >> logfile, "line number of inputfile:\t%s" % bamcount(bamfile)
    print >> logfile, "line number of inputfile after removing duplicated reads:\t%s" % bamcount(nodupbamfile)
    print >> logfile, "line number of inputfile after removing duplicated and multiple mapped reads:\t%s" % bamcount(nodupuniqbamfile)
    logfile.close()
    
    os.remove(bamfile)
    os.remove(nodupbamfile)
    
    endtime = time.time()
    print "it takes %.3f minutes or %.3f seconds" % ((endtime-starttime)/60, endtime -starttime)        