#author: zjuwhw
#date: 2014-12-25
USAGE='''extract_signal_from_bigwig.py --- extract signals from bigwig file accoarding to the bed file

USAGE:
    python %s [--type=#] [--strand=INT] [--nbins=#] [--range=#] [--output=<file>] bigwigfile bedfile(or stdin/-) 

#there are 4 types for this program
#type 1: the bed file is a region, output a average value appending echo line
#type 2: the bed file is a region and the nbins parameter is needed, output an aggregation plot which divide the region into nbins
#type 3: the bed file is a point and the range parameter is needed, output an aggregation plot which aggregate the signal in each position around the point
#type 4: the bed file is a point and the needed parameters are both nbins and range, output a heatmap plot
#--strand is just for type 2, 3, 4; and the INT is NUMBER of the strand column.
#input a bed formate file(the first three colnums are chr, str and end,respectively)

NOTE:
#the output value = real value/ row number of bedfile/ average tag density of bigWig file
#this program is based on the bigWigSummary tool, so the bigWigSummary must in the PATH
#if bedfile is stdin or -, the input is the stdin
#the bigwig file can be the /c/wanghw/annotation/conservation_score/genome.phastCons46way.bw when calculate conversion score.
#it is better to set the output parameter, because the bigWigSummary also output some contents into stdout
'''
import os,sys,getopt,time

def ossystemresult(command):
    fp=os.popen(command,"r")
    return fp.read()


def type1(bedfp,outfp,bigwigfile):
    for line in bedfp:
        linelist=line.rstrip().split("\t")
        cmd1="bigWigSummary %s %s %s %s %d" %(bigwigfile,linelist[0],linelist[1],linelist[2],1)
        score1=ossystemresult(cmd1)
        cmd2="bigWigSummary -type=coverage %s %s %s %s %d" %(bigwigfile,linelist[0],linelist[1],linelist[2],1)
        cover=ossystemresult(cmd2)
        if score1=="n/a" or score1=="" :
            score=0
        else:
            score=float(score1)*float(cover)
        print >> outfp, "%s\t%f" % (line.rstrip(),score)
        #print cmd1,cmd2,score1,cover,score
    bedfp.close()
    outfp.close()

def type2(bedfp,outfp,bigwigfile,nbins,nstrand):
    summary=[0]*nbins
    avr_tag_density=ossystemresult('''bigWigInfo %s |awk -F " " '$1=="mean:"{print $2}' ''' % bigwigfile)
    nlines=0
    for line in bedfp:
        nlines+=1
        linelist=line.rstrip().split("\t")
        cmd1="bigWigSummary %s %s %s %s %d" %(bigwigfile,linelist[0],linelist[1],linelist[2],nbins)
        score1=ossystemresult(cmd1)
        score1_list=score1.rstrip().split("\t")
        cmd2="bigWigSummary -type=coverage %s %s %s %s %d" %(bigwigfile,linelist[0],linelist[1],linelist[2],nbins)
        cover=ossystemresult(cmd2)        
        cover_list=cover.rstrip().split("\t")
        for i1 in range(nbins):
            if not nstrand:
                i2=i1
            else:
                if linelist[nstrand-1]=="+":
                    i2=i1
                elif linelist[nstrand-1]=="-":
                    i2=-(i1+1)
            try:
                score1_i=score1_list[i1]
                cover_i=cover_list[i1]
                summary[i2]+=float(score1_i)*float(cover_i)
            except:
                summary[i2]+=0
    for j in range(nbins):
        print >> outfp, "%d\t%.8f" % (j+1,float(summary[j])/float(nlines)/float(avr_tag_density))
    bedfp.close()
    outfp.close()
def type3(bedfp,outfp,bigwigfile,rangelength,nstrand):
    summary=[0]*(2*rangelength+1)
    avr_tag_density=ossystemresult('''bigWigInfo %s |awk -F " " '$1=="mean:"{print $2}' ''' % bigwigfile)
    nlines=0
    for line in bedfp:
        nlines+=1
        linelist=line.rstrip().split("\t")
        if int(linelist[2])-int(linelist[1]) !=1 :
            print "type 4 need input a point bed file!!"
            sys.exit(1)
        start=int(linelist[1])-rangelength
        end=int(linelist[2])+rangelength
        cmd1="bigWigSummary %s %s %d %d %d" %(bigwigfile,linelist[0],start,end,2*rangelength+1)
        score1=ossystemresult(cmd1)
        score1_list=score1.rstrip().split("\t")
        for i1 in range(2*rangelength+1):
            if not nstrand:
                i2=i1
            else:
                if linelist[nstrand-1]=="+":
                    i2=i1
                elif linelist[nstrand-1]=="-":
                    i2=-(i1+1)
            try:
                score1_i=score1_list[i1]
                summary[i2]+=float(score1_i)
            except:
                summary[i2]+=0
    #print nlines,summary,avr_tag_density,type(nlines),type(summary),len(summary),type(avr_tag_density),float(avr_tag_density)
    for j in range(2*rangelength+1):
        #print j-rangelength,float(summary[j])/float(nlines)/float(avr_tag_density)
        print >> outfp, "%d\t%.8f" % (j-rangelength,float(summary[j])/float(nlines)/float(avr_tag_density))
    bedfp.close()
    outfp.close()
def type4(bedfp,outfp,bigwigfile,nbins,rangelength,nstrand):
    avr_tag_density=ossystemresult('''bigWigInfo %s |awk -F " " '$1=="mean:"{print $2}' ''' % bigwigfile)
    for line in bedfp:
        linelist=line.rstrip().split("\t")
        if int(linelist[2])-int(linelist[1]) !=1 :
            print "type 4 need input a point bed file!!"
            sys.exit(1)
        start=int(linelist[1])-rangelength
        end=int(linelist[2])+rangelength
        cmd1="bigWigSummary %s %s %d %d %d" %(bigwigfile,linelist[0],start,end,nbins)
        score1=ossystemresult(cmd1)
        score1_list=score1.rstrip().split("\t")
        cmd2="bigWigSummary -type=coverage %s %s %s %s %d" %(bigwigfile,linelist[0],start,end,nbins)
        cover=ossystemresult(cmd2)        
        cover_list=cover.rstrip().split("\t")
        summary_line=[0]*nbins
        for i1 in range(nbins):
            if not nstrand:
                i2=i1
            else:
                if linelist[nstrand-1]=="+":
                    i2=i1
                elif linelist[nstrand-1]=="-":
                    i2=-(i1+1)
            try:
                score1_i=score1_list[i1]
                cover_i=cover_list[i1]
                summary_line[i2]+=float(score1_i)*float(cover_i)/float(avr_tag_density)
            except:
                summary_line[i2]+=0
                #print type(score1_i) , type(cover_i), score1_i, cover_i
        #print summary_line,avr_tag_density
        #print "\t".join(map(str,summary_line))
        print >> outfp,"\t".join(map(str,summary_line))
    bedfp.close()
    outfp.close()

    
def aggragation_plot(output):
    rfilename=output+".r"
    rfilefp = open(rfilename,"w")
    print >> rfilefp,'''
    file="%s"
    data=read.table(file)
    png(paste(file,".png",sep=""))
    plot(data[,1],data[,2],type="l",main=file)
    dev.off()''' % output
    rfilefp.close()

if __name__ == '__main__':
    starttime=time.time()
    if len(sys.argv) < 2:
        print USAGE % sys.argv[0]
        sys.exit(1)
    
    opts, args = getopt.getopt(sys.argv[1:], "", ["type=","strand=","nbins=","range=","output="])        
    
    if len(args)!=2:
        print "both the bigwigfile and the bedfile are needed\n" 
        sys.exit(1)
    
    # defaults
    ntype=1
    nbins=1
    nstrand=False
    rangelength=2000
    outfp=sys.stdout
    output="name"
    
    for o,a in opts:
        if o == '--type':
            ntype = int(a)
        elif o == '--strand':
            nstrand = int(a)
        elif o == '--nbins':
            nbins = int(a)
        elif o== '--range':
            rangelength = int(a)
        elif o == '--output':
            outfp = open(a, "w")
            output = a
    
    bigwigfile=args[0]
    if args[1]=="-" or args[1]=="stdin":
        bedfp=sys.stdin
    else:
        bedfp=open(args[1])
    
    if ntype==1:
        type1(bedfp,outfp,bigwigfile)
    elif ntype==2:
        type2(bedfp,outfp,bigwigfile,nbins,nstrand)
        #aggragation_plot(output)
    elif ntype==3:
        type3(bedfp,outfp,bigwigfile,rangelength,nstrand)
        #aggragation_plot(output)
    elif ntype==4:
        type4(bedfp,outfp,bigwigfile,nbins,rangelength,nstrand)
    else:
        print "the type parameter is needed anytime"
    endtime=time.time()
    print "it takes %d seconds to run this program!" % (endtime-starttime)
    
    
    
    
    
    
    
    
    
    
    
    
