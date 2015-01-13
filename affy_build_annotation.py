#author: zjuwhw
#date: 2014-12-27
USAGE = '''affy_build_annotation.py --- download annotation files for affy microarray probe id, using the R bioconductor annotation db package

    python %s affy_name("all" means all the affy microarray available)

Notes:
These affy microarrays are available: "hgu133a","hgu133a2","hgu133b","hgu133plus2","hgu219","hgu95a","hgu95av2","hgu95b","hgu95c","hgu95d","hgu95e","u133aaofav2"
The output annoation includes "GeneID", "GeneSymbol", "GeneAccnum", "GeneRefseqid", "GeneEnsemblid", "GeneName"
All annotation files are built by relative R .db packages, except the array "u133aaofav2", which was transformed from the GEO GPL file and deposits on my github.
'''

Rcode = '''
source("http://bioconductor.org/biocLite.R")
biocLite("%s.db")
library(%s.db)
id = keys(%s.db)
GeneId = unlist(as.list(%sENTREZID[id]))
GeneSymbol = unlist(as.list(%sSYMBOL[id]))
GeneAccnum = unlist(as.list(%sACCNUM[id]))
GeneRefseqid = sapply(as.list(%sREFSEQ[id]), function(x) x[1])
GeneEnsemblid = sapply(as.list(%sENSEMBL[id]), function(x) x[1])
GeneName = unlist(as.list(%sGENENAME[id]))
data = cbind(ProbeId = id, GeneId, GeneSymbol, GeneAccnum, GeneRefseqid, GeneEnsemblid, GeneName)
write.table(data, "%s_ann.txt", quote=F, sep="\t", row.names=F)
'''

Rcode2 = '''
install.packages("%s", repos=NULL)
install.packages("%s", repos=NULL)
'''

Rcode3 = '''
cat(is.element('%s', installed.packages()[,1]))
'''
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
                               "hgu95e":"http://mbni.org/customcdf/19.0.0/refseq.download/hgu95ehsrefseqcdf_19.0.0.tar.gz",
                               "u133aaofav2":"http://mbni.org/customcdf/19.0.0/refseq.download/u133aaofav2hsrefseqcdf_19.0.0.tar.gz"}
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
                               "hgu95e":"http://mbni.org/customcdf/19.0.0/ensg.download/hgu95ehsensgcdf_19.0.0.tar.gz",
                               "u133aaofav2":"http://mbni.org/customcdf/19.0.0/ensg.download/u133aaofav2hsensgcdf_19.0.0.tar.gz"}

import os,sys,getopt,time

def ossystemresult(command):
    fp=os.popen(command,"r")
    return fp.read()

def check_package(packagename):
    f = open("check_package.r","w")
    print >> f, Rcode3 % packagename
    f.close()
    cmd = "Rscript check_package.r"
    Result = ossystemresult(cmd)
    os.remove("check_package.r")
    return Result
   
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
            if check_package(sys.argv[1]+".db") == "TRUE":
                print "The annotation of affy microarray %s has been installed" % sys.argv[1]
            else:
                if sys.argv[1] != "u133aaofav2":
                    download_affy(sys.argv[1])
                else:
                    cmd = "wget --no-check-certificate https://github.com/zjuwhw/Pipeline/raw/master/supp_data/hgu133aaofav2_ann.txt.gz"
                    os.system(cmd)
                    cmd2 = "gunzip hgu133aaofav2_ann.txt.gz"
                    os.system(cmd2)
                    os.rename("hgu133aaofav2_ann.txt","u133aaofav2_ann.txt")
            
            customcdf_refseq_package_name = customcdfpackage_refseq_url[sys.argv[1]].split("/")[-1].split("_")[0]
            customcdf_ensg_package_name = customcdfpackage_ensg_url[sys.argv[1]].split("/")[-1].split("_")[0]
            if check_package(customcdf_ensg_package_name) == "TRUE" and check_package(customcdf_refseq_package_name) == "TRUE":
                print  "The customcdfpackage of %s has been installed" % sys.argv[1]
            else:
                download_customcdfpackage(sys.argv[1])
    
    endtime=time.time()
    print "it takes %d seconds or %d minutes or %d hours to run this program!" % (endtime-starttime, (endtime-starttime)/60, (endtime-starttime)/3600)
        
    
