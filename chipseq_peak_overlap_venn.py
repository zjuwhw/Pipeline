#author: zjuwhw
#date: 2014-12-25
usage='''peak_overlap_venn.py --- make the overlap venn plot of 2-5 peak bed file

USAGE:
    python %s [--output=<file>] _p_s.bedfiles

#the input file is 2-5 _p_s.bed peak bed file
#--output: the output file name. Default: overlap_venn_plot
#the bedtools must be installed in the PATH, and the VennDiagram package must be instaled in the R software.
'''

import itertools,os,sys,getopt
def countFilelines(file):
    count=-1
    for count,line in enumerate(open(file)):
        pass
    count+=1
    return str(count)
    
def ossystemresult(command):
    fp=os.popen(command,"r")
    return fp.read()
        
def overlap_result(files):
    result=[]
    for i in range(len(files)):
        combination_list=list(itertools.combinations(range(len(files)),i+1))
        for combination_id in combination_list:
            if len(combination_id)==1:
                result.append("area"+str(combination_id[0]+1)+"="+countFilelines(files[combination_id[0]]))
            elif len(combination_id)==2:
                cmd="bedtools intersect -u -a %s -b %s |wc -l" %(files[combination_id[0]],files[combination_id[1]])
                result.append("n"+str(combination_id[0]+1)+str(combination_id[1]+1)+"="+ossystemresult(cmd))
            else:
                cmd="bedtools intersect -u -a %s -b %s " %(files[combination_id[0]],files[combination_id[1]])
                prefix="n"+str(combination_id[0]+1)+str(combination_id[1]+1)
                for m in range(2,len(combination_id)):
                    cmd+="|bedtools intersect -u -a - -b %s" % files[combination_id[m]]
                    prefix+=str(combination_id[m]+1)
                cmd=cmd+"|wc -l"
                result.append(prefix+"="+ossystemresult(cmd))
    return result
    

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print usage % sys.argv[0]
        sys.exit(1)
    elif len(sys.argv) <3 or len(sys.argv) >6:
        print "please input 2-5 peak files"
    opts, args = getopt.getopt(sys.argv[1:], "", ["output="] )
    
    # defaults
    output = "overlap_venn_plot_"
    database = {2:"pairwise",3:"triple",4:"quad",5:"quintuple"}
    rscript = '''
library(VennDiagram)
venn.plot <- draw.%s.venn(
%s
,category = c(%s),
cat.cex = 0.5,
cex = 2
)
pdf(file=%s, width=8, height=8)
grid.draw(venn.plot);
dev.off();
''' 
    
    files = args
    filenames = []
    filenames2 = []
    for filename in files:
        filenames.append(os.path.basename(filename).replace(".bed","").replace("_p_s","").replace("_uniq_nodup",""))
        filenames2.append('"'+os.path.basename(filename).replace(".bed","").replace("_p_s","").replace("_uniq_nodup","")+'"')
    #print files,filenames
    output += "_".join(filenames)
    resultlist = overlap_result(files)
    print filenames,resultlist
    
    for o,a in opts:
        if o == '--output':
            output = a
    
    f = open(output+".r","w")
    if len(files) ==2 :
        if int(resultlist[0].replace("area1=",""))<int(resultlist[1].replace("area2=","")):
            print >> f, rscript % (database[len(files)], ",".join(resultlist).replace("n12","cross.area")+",inverted = T", ",".join([filenames2[1],filenames2[0]]),'"'+output+'.pdf"')
        else:
            print >> f, rscript % (database[len(files)], ",".join(resultlist).replace("n12","cross.area"), ",".join(filenames2),'"'+output+'.pdf"')
    else:
        print >> f, rscript % (database[len(files)], ",".join(resultlist), ",".join(filenames2),'"'+output+'.pdf"')
    f.close()
    os.system("Rscript "+output+".r")
    #os.remove("tmp_venn.r")
    os.remove("Rplots.pdf")
