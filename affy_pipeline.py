#author: zjuwhw
#date: 2014-12-26
USAGE='''affy_array_pipeline.py --- affy microarray pipeline using R to do rma, mas5.0 and/or not customCDF normalization
USAGE:

    python %s [--affyname=#] [--output=#] [--ann_affy_path=#] [--ann_refseq_bed12=#] cel_RAW_dictionary

Default:
--output is the basename of cel_RAW_dictionary
These affy microarrays are available: "hgu133a","hgu133a2","hgu133b","hgu133plus2","hgu219","hgu95a","hgu95av2","hgu95b","hgu95c","hgu95d","hgu95e","u133aaofav2"
--ann_affy_path is the path of affy annotation. The affy annoation can be built using another python program "affy_build_annotation.py". The default ann_affy_path is "/c/wanghw/annotation/affy/" and ends with "_ann.txt"
--ann_refseq_bed12 is the path of refseq bed12 annotation files, which could be downloaded from ucsc website. The default path is "/c/wanghw/annotation/refseq_hg19_07292013.bed"
'''
import os,sys,getopt,time

affys = ["hgu133a","hgu133a2","hgu133b","hgu133plus2","hgu219","hgu95a","hgu95av2","hgu95b","hgu95c","hgu95d","hgu95e","u133aaofav2"]
customcdfname_refseq = {"hgu133a":"HGU133A_Hs_REFSEQ",
                      "hgu133a2":"HGU133A2_Hs_REFSEQ",
                      "hgu133b":"HGU133B_Hs_REFSEQ",
                      "hgu133plus2":"HGU133Plus2_Hs_REFSEQ",
                      "hgu219":"HGU219_Hs_REFSEQ",
                      "hgu95a":"HGU95A_Hs_REFSEQ",
                      "hgu95av2":"HGU95Av2_Hs_REFSEQ",
                      "hgu95b":"HGU95B_Hs_REFSEQ",
                      "hgu95c":"HGU95C_Hs_REFSEQ",
                      "hgu95d":"HGU95D_Hs_REFSEQ",
                      "hgu95e":"HGU95E_Hs_REFSEQ",
                      "u133aaofav2":"U133AAofAv2_Hs_REFSEQ"}
customcdfname_ensg = {"hgu133a":"HGU133A_Hs_ENSG",
                      "hgu133a2":"HGU133A2_Hs_ENSG",
                      "hgu133b":"HGU133B_Hs_ENSG",
                      "hgu133plus2":"HGU133Plus2_Hs_ENSG",
                      "hgu219":"HGU219_Hs_ENSG",
                      "hgu95a":"HGU95A_Hs_ENSG",
                      "hgu95av2":"HGU95Av2_Hs_ENSG",
                      "hgu95b":"HGU95B_Hs_ENSG",
                      "hgu95c":"HGU95C_Hs_ENSG",
                      "hgu95d":"HGU95D_Hs_ENSG",
                      "hgu95e":"HGU95E_Hs_ENSG",
                      "u133aaofav2":"U133AAofAv2_Hs_ENSG"}


def cel2txt(cel_raw_dir, affy_type, name):
    print "Step1: it's begin to do normalization ..."
    Rcode = '''
dir="%s"
library(affy)
Data=ReadAffy(celfile.path=dir)
eRMA=rma(Data)
write.exprs(eRMA,file="%s_rma.txt")
eMAS = mas5(Data)
write.exprs(eMAS, file="%s_mas.txt")

#customcdf_refseq
Data=ReadAffy(celfile.path=dir, cdfname = "%s")
eRMA=rma(Data)
write.exprs(eRMA,file="%s_rma_customCDF_refseq.txt")
eMAS = mas5(Data)
write.exprs(eMAS, file="%s_mas_customCDF_refseq.txt")

#customcdf_ensg
#Data=ReadAffy(celfile.path=dir, cdfname = "%s")
#eRMA=rma(Data)
#write.exprs(eRMA,file="%s_rma_customCDF_ensg.txt")
#eMAS = mas5(Data)
#write.exprs(eMAS, file="%s_mas_customCDF_ensg.txt")
'''
    f = open("cel2txt_%s.r" % name, "w")
    print >> f, Rcode % (cel_raw_dir, name, name, customcdfname_refseq[affy_type], name, name, customcdfname_ensg[affy_type], name, name)
    f.close()
    cmd = "Rscript cel2txt_%s.r" % name
    os.system(cmd)
    os.remove("cel2txt_%s.r" % name)

def txt2anntxt(name, ann_refseq_bed12, ann_affy_path):
    print "Step2: it's begin to annotation the probe or refseq id using the gene symbol ..."
    cmd1 = r''' awk -F '\t' -v OFS='\t' 'BEGIN{while((getline<"%s")>0)l[$1]=$3}NR>1{$1=$1"|"l[$1]}{print $0}'  %s > %s ''' % (ann_affy_path, name + "_rma.txt", name + "_rma.ann.txt")
    cmd2 = r''' awk -F '\t' -v OFS='\t' 'BEGIN{while((getline<"%s")>0)l[$1]=$3}NR>1{$1=$1"|"l[$1]}{print $0}'  %s > %s ''' % (ann_affy_path, name + "_mas.txt", name + "_mas.ann.txt")
    cmd3 = r''' awk -F '\t' -v OFS='\t' 'BEGIN{while((getline<"%s")>0)l[$5]=$4}NR>1{split($1,x,".");if(x[1] in l){$1=$1"|"l[x[1]]}else{$1=$1"|NA"}}{print $0}'  %s > %s ''' % (ann_refseq_bed12, name + "_rma_customCDF_refseq.txt", name + "_rma_customCDF_refseq.ann.txt" )
    cmd4 = r''' awk -F '\t' -v OFS='\t' 'BEGIN{while((getline<"%s")>0)l[$5]=$4}NR>1{split($1,x,".");if(x[1] in l){$1=$1"|"l[x[1]]}else{$1=$1"|NA"}}{print $0}'  %s > %s ''' % (ann_refseq_bed12, name + "_mas_customCDF_refseq.txt", name + "_mas_customCDF_refseq.ann.txt" )
    #print cmd1
    #print cmd2
    #print cmd3
    #print cmd4
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    os.system(cmd4)

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
    
    opts, args = getopt.getopt(sys.argv[1:], "", ["affyname=","output=","ann_affy_path=","ann_refseq_bed12="])
    # defaults
    ann_refseq_bed12 = "/c/wanghw/annotation/refseq_hg19_07292013.bed"
    cel_raw_dir = args[0]
    nname = os.path.basename(cel_raw_dir)
    for o,a in opts:
        if o == '--affyname':
            naffyname = a
            ann_affy_path = "/c/wanghw/annotation/affy/%s_ann.txt" % naffyname
        elif o == '--output':
            nname = a
        elif o == '--ann_affy_path':
            ann_affy_path = a
        elif o == '--ann_refseq_bed12':
            ann_refseq_bed12 = a
            
    
    if naffyname not in affys:
        print "The %s is not in the list as below:" % naffyname
        print ", ".join(affys)
        sys.exit(1)
    
    cel2txt(cel_raw_dir, naffyname, nname)
    txt2anntxt(nname, ann_refseq_bed12, ann_affy_path)
    rowcollapse(nname + "_mas.ann.txt" )
    rowcollapse(nname + "_rma.ann.txt" )
    rowcollapse(nname + "_mas_customCDF_refseq.ann.txt" )
    rowcollapse(nname + "_rma_customCDF_refseq.ann.txt" )
    os.remove("%s_rma.txt" % nname)
    os.remove("%s_mas.txt" % nname)
    os.remove("%s_rma_customCDF_refseq.txt" % nname)
    os.remove("%s_mas_customCDF_refseq.txt" % nname)
    #os.remove("%s_rma_customCDF_ensg.txt" % nname)
    #os.remove("%s_mas_customCDF_ensg.txt" % nname)
    
    endtime=time.time()
    print "it takes %d seconds or %d minutes or %d hours to run this program!" % (endtime-starttime, (endtime-starttime)/60, (endtime-starttime)/3600)
    
