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

affys = ["hgu133a","hgu133a2","hgu133b","hgu133plus2","hgu219","hgu95a","hgu95av2","hgu95b","hgu95c","hgu95d","hgu95e"]

import os,sys,getopt,time

def download_affy(affy_name):
    print "The annotation of affy microarray %s is downloading ..." % affy_name
    f = open("affy_download_%s.r" % affy_name, "w")
    print >> f, Rcode % (affy_name, affy_name, affy_name, affy_name, affy_name, affy_name, affy_name, affy_name, affy_name, affy_name)
    f.close()
    cmd = "Rscript affy_download_%s.r"  % affy_name
    os.system(cmd)
    os.remove("affy_download_%s.r" % affy_name)

if __name__ == '__main__':
    starttime=time.time()
    if len(sys.argv) < 2:
        print USAGE % sys.argv[0]
        sys.exit(1)
    if sys.argv[1] == "all":
        for affy in affys:
            download_affy(affy)
    else:
        if sys.argv[1] not in affys:
            print "The %s is not in the list as below:" % sys.argv[1]
            print ", ".join(affys)
            sys.exit(1)
        else:
            download_affy(sys.argv[1])
    
    endtime=time.time()
    print "it takes %d seconds or %d minutes or %d hours to run this program!" % (endtime-starttime, (endtime-starttime)/60, (endtime-starttime)/3600)
        
    