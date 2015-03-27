#author: zjuwhw
#date: 2014-12-26
USAGE='''peak_calling.py --- peak calling for single-end ChIP-seq bam file using MACS and Hpeak software, and output the bigWig and peak region files
USAGE:
    python %s [--tools="macs","hpeak"] [--pvalue=#] [--output=#] [--chromsize=#]sample_bam_file control_bam_file(optional)

Default:
    This pipeline intergated two peak calling software, macs and hpeak. The tools parameter musted be set.
    For macs, the default p-value is 1e-8, and for hpeak, the default p-value is 1e-4
    For macs, the output are bigWig file and _p_s.bed file; for hpeak, the output is hpeak default output, further process should be made by yourself.
    the chromesize file is /c/wanghw/annotation/hg19.chrom.sizes, which could be downloaded from ucsc website.

Note:
    This pipeline needs the software macs, hpeak and ucsc's utility wigToBigWig installed in the PATH
    This pipeline was only tested using human data. You must pay more attention dealing with other species' data
'''

import os,sys,getopt,time

def macs_xls2psbed(macsxlsfile):
    fxls = open(macsxlsfile)
    fpsbed = open(macsxlsfile.replace("_peaks.xls","_p_s.bed"),"w")
    n = 0
    for line in fxls:
        if line.startswith("chr") and (not line.startswith("chr\t")):
            n = n+1
            linelist = line.rstrip().split("\t")
            print >> fpsbed, "%s\t%d\t%d\t%s\t%d" % (linelist[0], int(linelist[1])-1, int(linelist[2]), "MACS_peak_%s" % n, int(linelist[1])-1+int(linelist[4]) )
    fxls.close()
    fpsbed.close()
def macs_wig2bigWig(name, chromsizefile):
    cmd = "for gzwigfile in %s_MACS_wiggle/treat/*wig.gz;do gunzip -c $gzwigfile| awk 'NR>1{print $0}' >> %s;done" % (name, name + ".wig")
    os.system(cmd)
    cmd2 = "wigToBigWig -clip %s %s %s" % (name + ".wig", chromsizefile, name + ".bigWig")
    os.system(cmd2)
    cmd3 = "rm -rf %s_MACS_wiggle %s.wig" % (name, name)
    os.system(cmd3)
def macs_nocontrol(samplebamfile, pvalue, name, chromsize):
    cmd = "macs14 -t " + samplebamfile + " -p " + str(pvalue) + " -n " + name + " -w --space=10"
    os.system(cmd)
    cmd_rm="rm -f *r *negative*xls *peaks.bed *summits.bed *model.r"
    os.system(cmd_rm)
    macs_xls2psbed(name + "_peaks.xls")
    macs_wig2bigWig(name, chromsize)
def macs_havecontrol(samplebamfile, controlbamfile, pvalue, name, chromsize):
    cmd = "macs14 -t " + samplebamfile + " -c " + controlbamfile + " -p " + str(pvalue) + " -n " + name + " -w --space=10"
    os.system(cmd)
    cmd_rm="rm -f *r *negative*xls *peaks.bed *summits.bed *model.r"
    os.system(cmd_rm)
    macs_xls2psbed(name + "_peaks.xls")
    macs_wig2bigWig(name, chromsize)
def bam2hpeakbed(bamfile):
    cmd="bedtools bamtobed -i %s | cut -f1-3,6 > %s" %(bamfile, os.path.basename(bamfile).replace(".bam",".hpeak.bed"))
    os.system(cmd)
def hpeakout2psbed(name):
    cmd = '''awk '{if($1>=1 && $1<=22){$1="chr"$1}else if($1==23){$1="chrX"}else if($1==24){$1="chrY"};$2=int($2)-1;$5=int($5)+$2;$4="hpeak_peak_"NR;print $1"\t"$2"\t"$3"\t"$4"\t"$5}' %s >> %s ''' % (name+".hpeak.out", name+".hpeak_p_s.bed")
    #print cmd
    os.system(cmd)
def hpeak_nocontrol(samplebamfile, pvalue, name):
    bam2hpeakbed(samplebamfile)
    f1=open(name + "_treat.txt","w")
    print >> f1, os.path.basename(samplebamfile).replace(".bam",".hpeak.bed")
    f1.close()
    cmd="perl /c/wanghw/software/HPeak-2.1/HPeak.pl -format BED -t %s_treat.txt -n %s -fmin 100 -fmax 300 -w 25 -s %s" %(name, name, pvalue)
    os.system(cmd)
    os.remove(name + "_treat.txt")
    os.remove(os.path.basename(samplebamfile).replace(".bam",".hpeak.bed"))
    hpeakout2psbed(name)
def hpeak_havecontrol(samplebamfile, controlbamfile, pvalue, name):
    bam2hpeakbed(samplebamfile)
    bam2hpeakbed(controlbamfile)
    f1=open(name + "_treat.txt","w")
    f2=open(name + "_control.txt","w")
    print >> f1, os.path.basename(samplebamfile).replace(".bam",".hpeak.bed")
    print >> f2, os.path.basename(controlbamfile).replace(".bam",".hpeak.bed")
    f1.close()
    f2.close()
    cmd="perl /c/wanghw/software/HPeak-2.1/HPeak.pl -format BED -t %s_treat.txt -c %s_control.txt -n %s -fmin 100 -fmax 300 -w 25 -s %s" %(name, name, name, pvalue)
    os.system(cmd)
    os.remove(name + "_treat.txt")
    os.remove(name + "_control.txt")
    os.remove(os.path.basename(samplebamfile).replace(".bam",".hpeak.bed"))
    os.remove(os.path.basename(controlbamfile).replace(".bam",".hpeak.bed"))
    hpeakout2psbed(name+"_treated")

if __name__ == '__main__':
    starttime=time.time()
    if len(sys.argv) < 2:
        print USAGE % sys.argv[0]
        sys.exit(1)
    
    opts, args = getopt.getopt(sys.argv[1:], "", ["tools=","pvalue=","output="])
    if len(args) > 2:
        print "you can only input 1-2 bam file"
        sys.exit(1)
    if (('--tools','macs') not in opts) and (('--tools','hpeak') not in opts):
        print "the parameter tools must be set, and only the macs or hpeak could be used"
        sys.exit(1)
    
    # defaults
    def_pval={"macs": "1e-8", "hpeak": "1e-4"}
    nname = os.path.basename(args[0]).replace(".bam","")
    chrom_size = "/c/wanghw/annotation/hg19.chrom.sizes"
    
    for o,a in opts:
        if o == '--tools':
            ntype = a
            npvalue = def_pval[a]
        elif o == '--pvalue':
            npvalue = a
        elif o == '--output':
            nname = a
        elif o == '--chromsize':
            chrom_size = a
    print r'''#########''', args[0], npvalue, nname, chrom_size, r'''##############'''
    
    if len(args) == 2 and ntype == "macs":
        macs_havecontrol(args[0], args[1], npvalue, nname, chrom_size)
    elif len(args) == 1 and ntype == "macs":
        macs_nocontrol(args[0], npvalue, nname, chrom_size)
    elif len(args) == 2 and ntype == "hpeak":
        hpeak_havecontrol(args[0], args[1], npvalue, nname)
    elif len(args) == 1 and ntype == "hpeak":
        hpeak_nocontrol(args[0], npvalue, nname)    
    endtime=time.time()
    print "it takes %d seconds or %d minutes or %d hours to run this program!" % (endtime-starttime, (endtime-starttime)/60, (endtime-starttime)/3600)



