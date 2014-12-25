##author: zjuwhw
##time: 2014-12-25
import os,sys
import os.path
import time
def countFilelines(file):
    count=0
    f=open(file,"r")
    for line in f:
        count += 1
    return count

def sra2fastq(srafile):
    if os.path.exists(srafile.replace(".sra",".fastq"))  :
        print "The file "+srafile.replace(".sra",".fastq")+" has exist"
    elif os.path.exists(srafile.replace(".sra",".fq")):
        print "The file "+srafile.replace(".sra",".fq")+" has exist"
    else:
        cmd="fastq-dump "+srafile
        os.system(cmd)

def fastq2bam(fastqfile):
    postfix="."+os.path.basename(fastqfile).split(".")[-1]
    ref="/d/database/hg19/bwaindex/hg19.fa"
    if os.path.exists(fastqfile.replace(postfix,".bam"))  :
        print "The file "+fastqfile.replace(postfix,".bam")+" has exist"
    else:
        cmd1="bwa aln -t 4 "+ref+" "+fastqfile+" > "+fastqfile.replace(postfix,".sai")
        cmd2="bwa samse "+ref+" "+fastqfile.replace(postfix,".sai")+" "+fastqfile+" > "+fastqfile.replace(postfix,".sam")
        cmd3="samtools view -bS "+fastqfile.replace(postfix,".sam")+" > "+fastqfile.replace(postfix,"_unsorted.bam")
        cmd4="samtools sort "+fastqfile.replace(postfix,"_unsorted.bam")+" "+fastqfile.replace(postfix,"")
        cmd5="samtools index "+fastqfile.replace(postfix,"")+".bam"
        cmd6="rm -f *sam *sai *_unsorted.bam"
        os.system(cmd1)
        os.system(cmd2)
        os.system(cmd3)
        os.system(cmd4)
        os.system(cmd5)
        os.system(cmd6)

def bam2uniq_nodup_bam(bamfile):
    postfix="."+os.path.basename(bamfile).split(".")[-1]
    if os.path.exists(bamfile.replace(postfix,"_uniq_nodup.bam")):
        print "The file "+bamfile.replace(postfix,"_uniq_nodup.bam")+" has exist"
    else:
        cmd1="samtools view -bq 1 "+bamfile+" > "+bamfile.replace(postfix,"_uniq.bam")
        cmd2="samtools rmdup -s "+bamfile.replace(postfix,"_uniq.bam")+" "+bamfile.replace(postfix,"_uniq_nodup.bam")
        os.system(cmd1)
        os.system(cmd2)    
        
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'No action specified.'
        sys.exit()
    if sys.argv[1].startswith('--'):
        option = sys.argv[1][2:]
        if option=="help":
            print '''\
This program is used for ChIP-seq data preprocess and alignment by bwa software;
usage:
    python chipseq_pipeline.py [srafile or fastqfile or bamfile]
Notes:
    The program "bwa","samtools","sratoolkit" must be installed and avaliable in the PATH.
    The path of reference genome must be set and the bwa index must be bulit.
'''

        else :
            print 'Unknown option.'
        sys.exit()
    else:
        for filename in sys.argv[1:]:
            start=time.time()
            postfix="."+os.path.basename(filename).split(".")[-1]
            sra2fastq(filename.replace(postfix,".sra"))
            if os.path.exists(filename.replace(postfix,".fastq")):
                fastq2bam(filename.replace(postfix,".fastq"))
            else:
                fastq2bam(filename.replace(postfix,".fq"))
            bam2uniq_nodup_bam(filename.replace(postfix,".bam"))
            end=time.time()
            print "The file "+filename+" has done using "+str(end-start)+" seconds"
