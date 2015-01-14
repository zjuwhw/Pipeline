#author: zjuwhw
#date: 2015-01-13
USAGE='''affy_ExonOrGene_array_pipeline.py --- affy Exon Or Gene microarray pipeline using R package oligo to do rma normalization
USAGE:

    python %s [--affyname=#] [--output=#] [--ann_affy_path=#] cel_RAW_dictionary

Default:
--output is the basename of cel_RAW_dictionary
These affy microarrays are available: "huex10st","hugene10st","hugene11st","hugene20st","hugene21st"
--ann_affy_path is the path of affy annotation. The affy annoation can be built using another python program "affy_ExonOrGene_build_annotation.py".
The default ann_affy_path is "/c/wanghw/annotation/affy/" and ends with "_ann.txt"
Note:
Because the oligo package doesnot support the brainarray custom cdf of Exon/Gene ST arrays, this pipeline is only to do rma normalization with oligo package.
More details: https://groups.yahoo.com/neo/groups/customcdf/conversations/topics/294
'''
import os,sys,getopt,time

affys = ["huex10st","hugene10st","hugene11st","hugene20st","hugene21st"]
PDInfo ={"huex10st":"pd.huex.1.0.st.v2", "hugene10st":"pd.hugene.1.0.st.v1","hugene11st":"pd.hugene.1.1.st.v1", "hugene20st":"pd.hugene.2.0.st", "hugene21st":"pd.hugene.2.1.st"}

def cel2txt(cel_raw_dir, PDInfo_package_name, name):
    print "Step1: it's begin to do normalization ..."
    Rcode = '''
dir="%s"
library("oligo")
library("%s")
fns=list.files(path=dir,pattern='*.CEL.gz|*.CEL',full.names=T)
Data = read.celfiles(fns)
eRMA.core=rma(Data,target='core')
#featureData(eRMA.core)=getNetAffx(eRMA.core, 'transcript')
write.exprs(eRMA.core,file="%s_rma_core.txt")
#eRMA.probeset=rma(Data,target='probeset')
#featureData(eRMA.probeset)=getNetAffx(eRMA.probeset, 'probeset')
#write.exprs(eRMA.probeset,file="%s_rma_probeset.txt")
'''
    f = open("cel2txt_%s.r" % name, "w")
    print >> f, Rcode % (cel_raw_dir, PDInfo_package_name, name, name)
    f.close()
    cmd = "Rscript cel2txt_%s.r" % name
    os.system(cmd)
    os.remove("cel2txt_%s.r" % name)

def txt2anntxt(name, ann_affy_path):
    print "Step2: it's begin to annotation the probe or refseq id using the gene symbol ..."
    cmd1 = r''' awk -F '\t' -v OFS='\t' 'BEGIN{while((getline<"%s")>0)l[$1]=$3}NR>1{$1=$1"|"l[$1]}{print $0}'  %s > %s ''' % (ann_affy_path, name + "_rma_core.txt", name + "_rma_core.ann.txt")
    #print cmd1
    os.system(cmd1)

def rowcollapse(ann_file):
    print "Step3: it's begin to collapse the row with same gene symbol using the highest expression probe or refseq id ..."
    cmd1 = r'''head -1 %s > %s ''' % (ann_file, ann_file.replace("ann.txt","HighestUniqSymbols.txt"))
    cmd2 = r'''awk 'NR!=1{split($1,x,"|");$1=x[2];if($1!="NA"){n=0;for(i=2;i<=NF;i++){n=n+$i};print n,$0}}' %s | sort -k2,2 -k1,1nr|awk 'BEGIN{a="A"}$2!=a{printf $2;for(i=3;i<=NF;i++){printf "\t"$i};print "";a=$2}' >> %s '''% (ann_file, ann_file.replace("ann.txt","HighestUniqSymbols.txt"))
    #print cmd1
    #print cmd2
    os.system(cmd1)
    os.system(cmd2)

    
if __name__ == '__main__':
    starttime=time.time()
    if len(sys.argv) < 2:
        print USAGE % sys.argv[0]
        sys.exit(1)
    
    opts, args = getopt.getopt(sys.argv[1:], "", ["affyname=","output=","ann_affy_path="])
    # defaults
    #ann_refseq_bed12 = "/c/wanghw/annotation/refseq_hg19_07292013.bed"
    cel_raw_dir = args[0].rstrip("/")
    nname = os.path.basename(cel_raw_dir)
    for o,a in opts:
        if o == '--affyname':
            naffyname = a
            ann_affy_path = "/c/wanghw/annotation/affy/%stranscriptcluster_ann.txt" % naffyname
        elif o == '--output':
            nname = a
        elif o == '--ann_affy_path':
            ann_affy_path = a
    
    if naffyname not in affys:
        print "The %s is not in the list as below:" % naffyname
        print ", ".join(affys)
        sys.exit(1)
    
    cel2txt(cel_raw_dir, PDInfo[naffyname], nname)
    txt2anntxt(nname, ann_affy_path)
    rowcollapse(nname + "_rma_core.ann.txt" )
    os.remove("%s_rma_core.txt" % nname)
    #os.remove("%s_rma.txt" % nname)
    #os.remove("%s_mas.txt" % nname)
    #os.remove("%s_rma_customCDF_refseq.txt" % nname)
    #os.remove("%s_mas_customCDF_refseq.txt" % nname)
    #os.remove("%s_rma_customCDF_ensg.txt" % nname)
    #os.remove("%s_mas_customCDF_ensg.txt" % nname)
    
    endtime=time.time()
    print "it takes %d seconds or %d minutes or %d hours to run this program!" % (endtime-starttime, (endtime-starttime)/60, (endtime-starttime)/3600)
    
