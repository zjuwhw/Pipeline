USAGE='''fimo.py---used to run fimo programs in meme_suit for motif scaning
, caculate the numbers of motif sites in each peak
and output motif site position distribution centered around summit

USAGE:
    python %s [--output=#] [--RegionCenteredSummit=#] motiffile _p_s.bed option

#input: motif file and the _p_s.bed
#the motif file is in MEME format.It can be transfected to MEME format by meme-get-motif tool for the meme.txt file
#the p_s.bed file: the peak region is used to scan motif by FIMO and caculate the motif numbers in each peak; the summit region 
#defaults: 
#--output: the real name of motif file name and _p_s.bed file name
#--RangeAroundSummit: 200bp, which means the 201bp region +/- 100bp centered around summit
#option is the fimo option paramter, such as --thresh 1e-4

Note:
The tools, fimo and bedtools, are needed in $PATH
When cacualting the distance between motif site and summit, the motif is consided as a point in the middle.
'''

import os,sys,getopt,time,os.path,re,itertools, glob 

def ossystemresult(command):
    fp=os.popen(command,"r")
    return fp.read()
def getfasta_centor(psbed,d,name):
    cmd="awk '{print $1,$5-1,$5}' %s |bedtools slop -i - -g /c/wanghw/annotation/hg19.genome -b %d|bedtools getfasta -fi /d/database/hg19/hg19.fa -bed - -fo %s_tmp.fa" %(psbed, int(d/2), name)
    os.system(cmd)
def getfasta_region(psbed):
    cmd="bedtools getfasta -fi /d/database/hg19/hg19.fa -bed %s -fo %s_tmp.0-based.fa" % (psbed, name)
    os.system(cmd)
    cmd2=''' awk -F '\\t' -v OFS='\\t' '$1~/^>/{split($1,x,":|-");$1=x[1]":"(x[2]+1)"-"x[3]}{print $0}' %s_tmp.0-based.fa > %s_tmp.1-based.fa ''' % (name, name)
    os.system(cmd2)
    os.remove(name+"_tmp.0-based.fa")
def fimo(name,motiffile,option):
    cmd="fimo %s --parse-genomic-coord -oc %s %s %s_tmp.1-based.fa" % (option, name, motiffile, name)
    print cmd
    os.system(cmd)
    cmd2 = r'''awk  -F '\t' -v OFS='\t' 'NR==1{header="#"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}NR!=1{if(!($1 in container)){print header > "%s"$1"_fimo.bed";container[$1]=$1};print $2"\t"$3-1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9 >> "%s"$1"_fimo.bed"}' %s ''' % (name + "/", name + "/", name + "/" + "fimo.txt")
    os.system(cmd2)
    os.remove(name+"_tmp.1-based.fa")

def deletelaststring(dictory):
    return os.path.basename(dictory).rstrip("_fimo.bed")
def motifnumber(name,psbed):
    regex = name + "/*_fimo.bed"
    cmd = "cat %s" % psbed
    motifbedfiles = glob.glob(regex)
    for motifbedfile in motifbedfiles:
        cmd += "|bedtools intersect -c -a - -b %s " % (motifbedfile)
    cmd += r'''|awk 'BEGIN{print "chr\tstr\tend\tpeakid\tsummit\t%s"}{print $0}'> %s/motif_numeber_per_peak.txt''' % ("\\t".join(map(deletelaststring,motifbedfiles)), name)
    os.system(cmd)
    rscript = '''
data = read.table("%s/motif_numeber_per_peak.txt", header=T, sep="\t")
data = data[,6:ncol(data)]
library(ggplot2)
library(reshape2)
data.melt = melt(data)
data.melt[data.melt[,"value"] >= 6 ,"value"] = ">=6"
data.melt$value <- as.character(data.melt$value)
p <- ggplot(data.melt)+
  geom_histogram(aes(x=variable, fill=as.character(value)), position="fill")+
  xlab("motif")+
  ylab("percentage")+
  theme(panel.background=element_rect(colour="white",color="white",fill=NA),panel.border=element_rect(colour="black",fill=NA),legend.title=element_blank(),axis.text=element_text(colour="black"),panel.grid.minor = element_line(colour = NA),panel.grid.major = element_line(colour = NA),axis.text.x  = element_text(angle=30, vjust=0.5))
ggsave(filename="%s/motif_number_per_peak.pdf", plot=p)
''' % (name, name)
    f = open(name+"_rscript.r", "w")
    print >>f, rscript
    f.close()
    os.system("Rscript %s_rscript.r" % name)
    os.remove("%s_rscript.r" % name)

def allmotif2distance(name, psbed, d):
    regex = name + "/*_fimo.bed"
    motifbedfiles = glob.glob(regex)
    for motifbedfile in motifbedfiles:
        cmd = r'''bedtools intersect -wao -a %s -b %s |awk '{print "%s\t"int(($2+$3)/2)-$(NF-1)}' >> %s/allmotif2distance.txt ''' % (motifbedfile, psbed, deletelaststring(motifbedfile), name)
        os.system(cmd)
    rscript = '''
data = read.table("%s/allmotif2distance.txt", sep="\t")
colnames(data) = c("motif","motif2summit")
library(ggplot2)
p <- ggplot(data)+
  geom_density(aes(x = motif2summit, colour=motif))+
  xlim(-%d,%d)+
  ggtitle("Motif position distribution")
ggsave(filename="%s/allmotif2distance.pdf", plot=p)
''' % (name, int(d/2), int(d/2), name)
    f = open(name+"_rscript.r", "w")
    print >> f, rscript
    f.close()
    os.system("Rscript %s_rscript.r" % name)
    os.remove("%s_rscript.r" % name)
    os.remove(name + "/allmotif2distance.txt")


if __name__ == '__main__':
    
    starttime=time.time()
    if len(sys.argv) < 2:
        print USAGE % sys.argv[0]
        sys.exit(1)
        
    opts, args = getopt.getopt(sys.argv[1:], "", ["output=","RegionCenteredSummit="])
    
    # defaults
    option="--thresh 1e-4"
    motiffile=args[0]
    psbed = args[1]
    name="fimo_output_"+os.path.basename(motiffile)+"_"+os.path.basename(psbed)
    d=200
    
    for o,a in opts:
        if o == '--output':
            name = "fimo_output_"+a
        elif o == '--RegionCenteredSummit':
            d = int(a)
    
    if len(args)!= 2:
        option = " ".join(args[2:])
    
    getfasta_region(psbed)
    fimo(name,motiffile,option)
    motifnumber(name,psbed)
    allmotif2distance(name, psbed, d)
    endtime = time.time()
    print "it takes %.3f secondes or %.3f minutes" % (endtime -starttime, (endtime-starttime)/60)
