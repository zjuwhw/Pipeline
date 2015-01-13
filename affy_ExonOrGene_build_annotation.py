#author: zjuwhw
#date: 2015-01-13
USAGE = '''affy_ExonOrGene_build_annotation.py --- download annotation files for affy Exon Or Gene microarray transcriptcluster id, using the R bioconductor annotation db package

    python %s affy_name("all" means all the affy microarray available)

Notes:
These affy Exon Or Gene microarrays are available: "huex10st","hugene10st","hugene11st","hugene20st","hugene21st"
The output transcriptcluster annoation includes "GeneID", "GeneSymbol", "GeneAccnum", "GeneRefseqid", "GeneEnsemblid", "GeneName"
All annotation files are built by relative R transcriptcluster.db packages
'''

Rcode_check_package = '''
cat(is.element('%s', installed.packages()[,1]))
'''
Rcode_install_biocLite = '''
source("http://bioconductor.org/biocLite.R")
biocLite("%s")
'''
Rcode_install_local = '''
install.packages("%s", repos=NULL)
'''
Rcode_ann_transcriptcluster = '''
library(%stranscriptcluster.db)
id = keys(%stranscriptcluster.db)
GeneId = unlist(as.list(%stranscriptclusterENTREZID[id]))
GeneSymbol = unlist(as.list(%stranscriptclusterSYMBOL[id]))
GeneAccnum = sapply(as.list(%stranscriptclusterACCNUM[id]), function(x) x[1])
GeneRefseqid = sapply(as.list(%stranscriptclusterREFSEQ[id]), function(x) x[1])
GeneEnsemblid = sapply(as.list(%stranscriptclusterENSEMBL[id]), function(x) x[1])
GeneName = unlist(as.list(%stranscriptclusterGENENAME[id]))
data = cbind(TranscriptclusterId = id, GeneId, GeneSymbol, GeneAccnum, GeneRefseqid, GeneEnsemblid, GeneName)
write.table(data, "%stranscriptcluster_ann.txt", quote=F, sep="\t", row.names=F)
'''

affys = ["huex10st","hugene10st","hugene11st","hugene20st","hugene21st"]
PDInfo ={"huex10st":"pd.huex.1.0.st.v2", "hugene10st":"pd.hugene.1.0.st.v1","hugene11st":"pd.hugene.1.1.st.v1", "hugene20st":"pd.hugene.2.0.st", "hugene21st":"pd.hugene.2.1.st"}
customcdfname_refseq = {"huex10st":"huex10st_Hs_REFSEQ",
                        "hugene10st":"hugene10st_Hs_REFSEQ",
                        "hugene11st":"hugene11st_Hs_REFSEQ",
                        "hugene20st":"hugene20st_Hs_REFSEQ",
                        "hugene21st":"hugene21st_Hs_REFSEQ"}
customcdfname_ensg = {"huex10st":"huex10st_Hs_ENSG",
                      "hugene10st":"hugene10st_Hs_ENSG",
                      "hugene11st":"hugene11st_Hs_ENSG",
                      "hugene20st":"hugene20st_Hs_ENSG",
                      "hugene21st":"hugene21st_Hs_ENSG"}
customcdfpackage_refseq_url = {"huex10st":"http://mbni.org/customcdf/19.0.0/refseq.download/huex10sthsrefseqcdf_19.0.0.tar.gz",
                               "hugene10st":"http://mbni.org/customcdf/19.0.0/refseq.download/hugene10sthsrefseqcdf_19.0.0.tar.gz",
                               "hugene11st":"http://mbni.org/customcdf/19.0.0/refseq.download/hugene11sthsrefseqcdf_19.0.0.tar.gz",
                               "hugene20st":"http://mbni.org/customcdf/19.0.0/refseq.download/hugene20sthsrefseqcdf_19.0.0.tar.gz",
                               "hugene21st":"http://mbni.org/customcdf/19.0.0/refseq.download/hugene21sthsrefseqcdf_19.0.0.tar.gz"}
customcdfpackage_ensg_url = {"huex10st":"http://mbni.org/customcdf/19.0.0/ensg.download/huex10sthsensgcdf_19.0.0.tar.gz",
                               "hugene10st":"http://mbni.org/customcdf/19.0.0/ensg.download/hugene10sthsensgcdf_19.0.0.tar.gz",
                               "hugene11st":"http://mbni.org/customcdf/19.0.0/ensg.download/hugene11sthsensgcdf_19.0.0.tar.gz",
                               "hugene20st":"http://mbni.org/customcdf/19.0.0/ensg.download/hugene20sthsensgcdf_19.0.0.tar.gz",
                               "hugene21st":"http://mbni.org/customcdf/19.0.0/ensg.download/hugene21sthsensgcdf_19.0.0.tar.gz"}

import os,sys,getopt,time

def ossystemresult(command):
    fp=os.popen(command,"r")
    return fp.read()
def check_package(packagename):
    f = open("check_package.r","w")
    print >> f, Rcode_check_package % packagename
    f.close()
    cmd = "Rscript check_package.r"
    Result = ossystemresult(cmd)
    os.remove("check_package.r")
    if Result == "TRUE":
        print "The package %s has been installed." % packagename
        return True
    elif Result == "FALSE":
        print "The package %s has not been installed." % packagename
        return False
def install_biocLite(packagename):
    print "It's going to install the package %s ..." % packagename
    f = open("install_biocLite.r","w")
    print >> f, Rcode_install_biocLite % packagename
    f.close()
    cmd = "Rscript install_biocLite.r"
    os.system(cmd)
    os.remove("install_biocLite.r")
def install_local(packagefilename):
    print "It's going to install the package %s ..." % packagefilename
    f = open("install_local.r","w")
    print >> f, Rcode_install_local % packagefilename
    f.close()
    cmd = "Rscript install_local.r"
    os.system(cmd)
    os.remove("install_local.r")
def install_customcdf(affy_name):
    cmd1 = "wget %s" % customcdfpackage_refseq_url[affy_name]
    cmd2 = "wget %s" % customcdfpackage_ensg_url[affy_name]
    os.system(cmd1)
    os.system(cmd2)
    install_local(customcdfpackage_refseq_url[affy_name].split("/")[-1])
    install_local(customcdfpackage_ensg_url[affy_name].split("/")[-1])
    os.remove(customcdfpackage_refseq_url[affy_name].split("/")[-1])
    os.remove(customcdfpackage_ensg_url[affy_name].split("/")[-1])
    
def ann_transcirptcluster(affy_name):
    f = open("ann_transcriptcluster_%s.r" % affy_name, "w")
    print >> f, Rcode_ann_transcriptcluster % (affy_name, affy_name, affy_name, affy_name, affy_name, affy_name, affy_name, affy_name, affy_name)
    f.close()
    cmd = "Rscript ann_transcriptcluster_%s.r"  % affy_name
    os.system(cmd)
    os.remove("ann_transcriptcluster_%s.r" % affy_name)    


if __name__ == '__main__':
    starttime=time.time()
    if len(sys.argv) < 2:
        print USAGE % sys.argv[0]
        sys.exit(1)
    if sys.argv[1] == "all":
        for affy in affys:
            PDInfo_package_name = PDInfo[affy]
            transcriptcluster_db_package_name = affy + "transcriptcluster.db"
            customcdf_refseq_package_name = customcdfpackage_refseq_url[affy].split("/")[-1].split("_")[0]
            customcdf_ensg_package_name = customcdfpackage_ensg_url[affy].split("/")[-1].split("_")[0]
            if not check_package(transcriptcluster_db_package_name):
                install_biocLite(transcriptcluster_db_package_name)
            if not check_package(PDInfo_package_name):
                install_biocLite(PDInfo_package_name)
            if (not check_package(customcdf_ensg_package_name)) or (not check_package(customcdf_refseq_package_name)):
                install_customcdf(affy)            
            ann_transcirptcluster(affy)
    else:
        affy = sys.argv[1]
        if affy not in affys:
            print "The %s is not in the list as below:" % sys.argv[1]
            print ", ".join(affys)
            sys.exit(1)
        else:
            PDInfo_package_name = PDInfo[affy]
            transcriptcluster_db_package_name = affy + "transcriptcluster.db"
            customcdf_refseq_package_name = customcdfpackage_refseq_url[affy].split("/")[-1].split("_")[0]
            customcdf_ensg_package_name = customcdfpackage_ensg_url[affy].split("/")[-1].split("_")[0]
            if not check_package(transcriptcluster_db_package_name):
                install_biocLite(transcriptcluster_db_package_name)
            if not check_package(PDInfo_package_name):
                install_biocLite(PDInfo_package_name)
            if (not check_package(customcdf_ensg_package_name)) or (not check_package(customcdf_refseq_package_name)):
                install_customcdf(affy)            
            ann_transcirptcluster(affy)
    endtime=time.time()
    print "it takes %d seconds or %d minutes or %d hours to run this program!" % (endtime-starttime, (endtime-starttime)/60, (endtime-starttime)/3600)
        
    
