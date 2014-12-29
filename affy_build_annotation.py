#author: zjuwhw
#date: 2014-12-27
USAGE = '''affy_build_annotation.py --- download annotation files for affy microarray probe id, using the R bioconductor annotation db package

    python %s affy_name("all" means all the affy microarray available)

Notes:
These affy microarrays are available: "hgu133a","hgu133a2","hgu133b","hgu133plus2","hgu219","hgu95a","hgu95av2","hgu95b","hgu95c","hgu95d","hgu95e"
The output annoation includes "GeneID", "GeneSymbol", "GeneAccnum", "GeneRefseqid", "GeneEnsemblid", "GeneName"
'''

Rcode = '''
source("http://bioconductor.org/biocLite.R")
biocLite("%s.db")
library(%s.db)
id = keys(%sCHR)
GeneId = unlist(as.list(%sENTREZID[id]))
GeneSymbol = unlist(as.list(%sSYMBOL[id]))
GeneAccnum = unlist(as.list(%sACCNUM[id]))
GeneRefseqid = sapply(as.list(%sREFSEQ[id]), function(x) x[1])
GeneEnsemblid = sapply(as.list(%sENSEMBL[id]), function(x) x[1])
GeneName = unlist(as.list(%sGENENAME[id]))
data = cbind(GeneId, GeneSymbol, GeneAccnum, GeneRefseqid, GeneEnsemblid, GeneName)
data = cbind(ProbeId = rownames(data), data)
write.table(data, "%s_ann.txt", quote=F, sep="\t", row.names=F)
'''

Rcode2 = '''
install.packages("%s", repos=NULL)
install.packages("%s", repos=NULL)
'''

affys = ["hgu133a","hgu133a2","hgu133b","hgu133plus2","hgu219","hgu95a","hgu95av2","hgu95b","hgu95c","hgu95d","hgu95e"]
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
                      "hgu95e":"HGU95E_Hs_REFSEQ"}
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
                      "hgu95e":"HGU95E_Hs_ENSG"}
customcdfpackage_refseq_url = {"hgu133a":"http://mbni.org/customcdf/19.0.0/refseq.download/hgu133ahsrefseqcdf_19.0.0.tar.gz",
                               "hgu133a2":"http://mbni.org/customcdf/19.0.0/refseq.download/hgu133a2hsrefseqcdf_19.0.0.tar.gz",
                               "hgu133b":"http://mbni.org/customcdf/19.0.0/refseq.download/hgu133bhsrefseqcdf_19.0.0.tar.gz",
                               "hgu133plus2":"http://mbni.org/customcdf/19.0.0/refseq.download/hgu133plus2hsrefseqcdf_19.0.0.tar.gz",
                               "hgu219":"http://mbni.org/customcdf/19.0.0/refseq.download/hgu219hsrefseqcdf_19.0.0.tar.gz",
                               "hgu95a":"http://mbni.org/customcdf/19.0.0/refseq.download/hgu95ahsrefseqcdf_19.0.0.tar.gz",
                               "hgu95av2":"http://mbni.org/customcdf/19.0.0/refseq.download/hgu95av2hsrefseqcdf_19.0.0.tar.gz",
                               "hgu95b":"http://mbni.org/customcdf/19.0.0/refseq.download/hgu95bhsrefseqcdf_19.0.0.tar.gz",
                               "hgu95c":"http://mbni.org/customcdf/19.0.0/refseq.download/hgu95chsrefseqcdf_19.0.0.tar.gz",
                               "hgu95d":"http://mbni.org/customcdf/19.0.0/refseq.download/hgu95dhsrefseqcdf_19.0.0.tar.gz",
                               "hgu95e":"http://mbni.org/customcdf/19.0.0/refseq.download/hgu95ehsrefseqcdf_19.0.0.tar.gz"}
customcdfpackage_ensg_url = {"hgu133a":"http://mbni.org/customcdf/19.0.0/ensg.download/hgu133ahsensgcdf_19.0.0.tar.gz",
                               "hgu133a2":"http://mbni.org/customcdf/19.0.0/ensg.download/hgu133a2hsensgcdf_19.0.0.tar.gz",
                               "hgu133b":"http://mbni.org/customcdf/19.0.0/ensg.download/hgu133bhsensgcdf_19.0.0.tar.gz",
                               "hgu133plus2":"http://mbni.org/customcdf/19.0.0/ensg.download/hgu133plus2hsensgcdf_19.0.0.tar.gz",
                               "hgu219":"http://mbni.org/customcdf/19.0.0/ensg.download/hgu219hsensgcdf_19.0.0.tar.gz",
                               "hgu95a":"http://mbni.org/customcdf/19.0.0/ensg.download/hgu95ahsensgcdf_19.0.0.tar.gz",
                               "hgu95av2":"http://mbni.org/customcdf/19.0.0/ensg.download/hgu95av2hsensgcdf_19.0.0.tar.gz",
                               "hgu95b":"http://mbni.org/customcdf/19.0.0/ensg.download/hgu95bhsensgcdf_19.0.0.tar.gz",
                               "hgu95c":"http://mbni.org/customcdf/19.0.0/ensg.download/hgu95chsensgcdf_19.0.0.tar.gz",
                               "hgu95d":"http://mbni.org/customcdf/19.0.0/ensg.download/hgu95dhsensgcdf_19.0.0.tar.gz",
                               "hgu95e":"http://mbni.org/customcdf/19.0.0/ensg.download/hgu95ehsensgcdf_19.0.0.tar.gz"}

import os,sys,getopt,time

def download_affy(affy_name):
    print "The annotation of affy microarray %s is downloading ..." % affy_name
    f = open("affy_download_%s.r" % affy_name, "w")
    print >> f, Rcode % (affy_name, affy_name, affy_name, affy_name, affy_name, affy_name, affy_name, affy_name, affy_name, affy_name)
    f.close()
    cmd = "Rscript affy_download_%s.r"  % affy_name
    os.system(cmd)
    os.remove("affy_download_%s.r" % affy_name)
    print "It is successful to download the nanotation of affy microarray %s" % affy_name
def download_customcdfpackage(affy_name):
    print "The customcdfpackage of %s is downloading ..." % affy_name
    cmd1 = "wget %s" % customcdfpackage_refseq_url[affy_name]
    cmd2 = "wget %s" % customcdfpackage_ensg_url[affy_name]
    os.system(cmd1)
    os.system(cmd2)
    f = open("customcdf_download_%s.r" % affy_name, "w")
    print >> f, Rcode2 % (customcdfpackage_refseq_url[affy_name].split("/")[-1], customcdfpackage_ensg_url[affy_name].split("/")[-1])
    f.close()
    cmd = "Rscript customcdf_download_%s.r" % affy_name
    os.system(cmd)
    os.remove("customcdf_download_%s.r" % affy_name)
    os.remove(customcdfpackage_refseq_url[affy_name].split("/")[-1])
    os.remove(customcdfpackage_ensg_url[affy_name].split("/")[-1])
    print "Ti is successful to install customcdfpackage of %s"  % affy_name
    


if __name__ == '__main__':
    starttime=time.time()
    if len(sys.argv) < 2:
        print USAGE % sys.argv[0]
        sys.exit(1)
    if sys.argv[1] == "all":
        for affy in affys:
            download_affy(affy)
            download_customcdfpackage(affy)
    else:
        if sys.argv[1] not in affys:
            print "The %s is not in the list as below:" % sys.argv[1]
            print ", ".join(affys)
            sys.exit(1)
        else:
            download_affy(sys.argv[1])
            download_customcdfpackage(sys.argv[1])
    
    endtime=time.time()
    print "it takes %d seconds or %d minutes or %d hours to run this program!" % (endtime-starttime, (endtime-starttime)/60, (endtime-starttime)/3600)
        
    