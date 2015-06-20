USAGE='''ngs_QCandTrimm.py--- Qualtiy Control and Trimming for single/pair-end ngs fastq files.
Quality Control for fastq before and after trimming using FastQC,
Qualtiy and adapter trimming using trimmomatic.

USAGE:
    python %s [--type=#] [--trimmomatic_option=#][--adapter=#] (read_1.fastq read_2.fastq / read.fastq / read.sra)

#input: single/pair-end fastq or sra files, it must obey the regular expression "(((_[12])?\.(fastq|fq)(\.gz)?)|(\.sra))$"
#output: a folder named ${prefix}_QCandTrim, including fastqc.html files and trimmed.fastq files (prefix: basename of input file(s))

#defaults:
--type: "SE" or "PE", the option must be set.
--trimmomatic_option: ILLUMINACLIP:$adapter:2:30:10 SLIDINGWINDOW:4:15 TRAILING:10 MINLEN:20, which means
1)Remove Illumina adapters provided in the adapter fasta file. Initially Trimmomatic will look for seed matches (16 bases) allowing maximally 2 mismatches. These seeds will be extended and clipped if in the case of paired end reads a score of 30 is reached (about 50 bases), or in the case of single end reads a score of 10, (about 17 bases)
2)Remove trailing low quality or N bases (below quality 10)
3)scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15.
4)Drop reads which are less than 20 bases long after these steps.
--adapter: the path of adapter fasta file, trimmomatic provides "TruSeq2-PE.fa" for GAII machines and "TruSeq3-PE.fa" or "TruSeq3-PE-2.fa" for HiSeq and MiSeq machines. Default is "/c/wanghw/software/rnaseq/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa".

Note:
The tools, fastqc and trimmomatic, sratoolkit are needed in $PATH
'''

import os, getopt, sys, time, os.path, re

def sra2fastq_SE(srafile, prefix):
    cmd = "fastq-dump -Z %s > %s" % (srafile, prefix+".fastq")
    os.system(cmd)
    return prefix+".fastq"
def sra2fastq_PE(srafile, prefix):
    cmd = "fastq-dump --split-files %s" % srafile
    os.system(cmd)
    return [prefix+"_1.fastq", prefix+"_2.fastq"]
def fastqc(fastqfile):
    cmd = "fastqc -o ./ %s" % (fastqfile)
    os.system(cmd)
    print "########fastqc command:\t" + cmd
def trim_PE(read1, read2, prefix, postfix, option):
    read1P_name=prefix+"_trimmedP_1"+postfix
    read1U_name=prefix+"_trimmedU_1"+postfix
    read2P_name=prefix+"_trimmedP_2"+postfix
    read2U_name=prefix+"_trimmedU_2"+postfix
    cmd = "trimmomatic PE %s %s %s %s %s %s %s" % (read1, read2, read1P_name, read1U_name, read2P_name, read2U_name, option)
    print "########Trimmomatic PE command:\t" + cmd
    os.system(cmd)
    return [read1P_name, read1U_name, read2P_name, read2U_name]
def trim_SE(read, prefix, postfix, option):
    readP_name = prefix + "_trimmed" + postfix
    cmd = "trimmomatic SE %s %s %s " % (read, readP_name, option)
    os.system(cmd)
    print "########Trimmomatic SE command:\t" + cmd
    return readP_name

if __name__ == '__main__':
    starttime=time.time()
    if len(sys.argv) < 2:
        print USAGE % sys.argv[0]
        sys.exit(1)
        
    # defaults    
    adapter="/c/wanghw/software/rnaseq/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa"
    option="ILLUMINACLIP:%s:2:30:10 SLIDINGWINDOW:4:15 TRAILING:10 MINLEN:20" % adapter
    
    opts, args = getopt.getopt(sys.argv[1:], "", ["type=","trimmomatic_option=","adapter="])        
    for o,a in opts:
        if o == "--type":
            datatype = a
        if o == "--adapter":
            adapter = a
            option = "ILLUMINACLIP:%s:2:30:10 SLIDINGWINDOW:4:15 TRAILING:10 MINLEN:20" % a
        if o == "--trimmomatic_option":
            option = a
    
    regular = "(((_[12])?\.(fastq|fq)(\.gz)?)|(\.sra))$"
    regular2 = "((\.(fastq|fq)(\.gz)?)|(\.sra))$"
    if not re.search(regular, args[0]):
        print '''the input file must obey the regular expression "(((_[12])?\.(fastq|fq)(\.gz)?)|(\.sra))$" '''
        sys.exit(1)
    else:
        prefix = re.sub(regular, "", os.path.basename(args[0]))
        postfix = re.search(regular2, os.path.basename(args[0])).group(0)    
    
    
    if not os.path.exists(prefix+"_QCandTrim"):
        os.mkdir(prefix+"_QCandTrim")
    if datatype == "SE":
        if postfix == ".sra":
            read = sra2fastq_SE(args[0], prefix)
            postfix = ".fastq"
        else:
            read = args[0]        
        fastqc(read)
        trimmedread = trim_SE(read, prefix, postfix, option)
        fastqc(trimmedread)
        os.system("mv %s*fastqc.zip %s*fastqc.html %s %s" % (prefix, prefix, trimmedread, prefix+"_QCandTrim"))
        
    elif datatype == "PE":
        if postfix == ".sra":
            read1, read2 = sra2fastq_PE(args[0], prefix)
            postfix = ".fastq"
        else:
            read1 = args[0]
            read2 = args[1]
        fastqc(read1)
        fastqc(read2)
        trimmedread1P, trimmedread1U, trimmedread2P, trimmedread2U = trim_PE(read1, read2, prefix, postfix, option)
        fastqc(trimmedread1P)
        fastqc(trimmedread2P)
        os.system("mv %s*fastqc.zip %s*fastqc.html %s %s %s %s %s" % (prefix, prefix, trimmedread1P, trimmedread1U, trimmedread2P, trimmedread2U, prefix+"_QCandTrim"))

    endtime = time.time()
    print "it takes %.3f secondes or %.3f minutes" % (endtime -starttime, (endtime-starttime)/60)
