#Author: zjuwhw
#date: 2014-12-25
USAGE='''motif_discovery.py---used for motif discovery of ChIP-seq datas, +/-100bp region around summit using meme-chip, homer or amd software.

USAGE:
    python %s [--tools="meme-chip","homer","amd"] [--sequence=#] [--motif_db=#] [--output=#] [--d=#] _p_s.bed optional(specific for the tools)
    
#input is the _p_s.bed
#defaults: 
 #--motif_db /c/wanghw/motif_database/Transfect_9.2.meme 
 #--output the real name of _p_s.bed file name
 #--sequence /d/database/hg19/hg19.fa.masked
 #--d 200bp, when d is "fulllength" means using amd and peak full length
 #optional:
  For meme-chip: -meme-nmotifs 5 -meme-minw 6 -meme-maxw 20 
  For homer: -mask -S 20 -len 8,10,12,14 -p 4 -size 200
  For amd: -T 50 -CO 0.6 -FC 1.2
Note:
The tools, meme-chip, amd, homer and bedtools, are needed in $PATH
'''
import os,sys,getopt,time

def ossystemresult(command):
    fp=os.popen(command,"r")
    return fp.read()
def getbed_centor(psfp, d):
	psfp.seek(0)
	tmpbedfp=open("tmp.bed","w")
	for line in psfp:
		linelist=line.rstrip().split("\t")
		print >> tmpbedfp, "%s\t%s\t%s\t%s" % (linelist[0],str(int(linelist[4])-1),linelist[4],linelist[3])
	tmpbedfp.close()
	cmd = ''' bedtools slop -i tmp.bed -g /c/wanghw/annotation/hg19.genome -b %d | awk 'BEGIN{OFS="\\t"}{print $0,".","+"; print $0, ".", "-"}' > tmp_homer.bed''' %(int(d/2))
	os.system(cmd)
	os.remove("tmp.bed")
def getfasta_centor(psfp, d, sequence):
	psfp.seek(0)
	tmpbedfp=open("tmp.bed","w")
	for line in psfp:
		linelist=line.rstrip().split("\t")
		print >>tmpbedfp,"%s\t%s\t%s\t%s" % (linelist[0],str(int(linelist[4])-1),linelist[4],linelist[3])
	tmpbedfp.close()
	cmd="bedtools slop -i tmp.bed -g /c/wanghw/annotation/hg19.genome -b %d|bedtools getfasta -fi %s -bed - -fo tmp.fa" %(int(d/2), sequence)
	os.system(cmd)
	os.remove("tmp.bed")
def getfasta_fulllength(psfp, sequence):
	psfp.seek(0)
	tmpbedfp=open("tmp.bed","w")
	for line in psfp:
		linelist=line.rstrip().split("\t")
		print >>tmpbedfp,"%s\t%s\t%s\t%s" % (linelist[0],linelist[1],linelist[2],linelist[3])
	tmpbedfp.close()
	cmd="bedtools getfasta -fi %s -bed %s -fo tmp.fa" %(sequence, "tmp.bed")
	os.system(cmd)
	os.remove("tmp.bed")
def getfasta_region(psfp):
	psfp.seek(0)
	tmpbedfp=open("tmp.bed","w")
	for line in psfp:
		linelist=line.rstrip().split("\t")
		print >> tmpbedfp,"%s\t%s\t%s\t%s" % (linelist[0],linelist[1],linelist[2],linelist[3])
	cmd="bedtools getfasta -fi %s -bed tmp.bed -fo tmp.fa" % (sequence)
	os.system(cmd)
	os.remove("tmp.bed")
def meme_chip(name, motif_db, optional):
	cmd = "meme-chip %s -o %s -db %s tmp.fa" % (optional, name, motif_db)
	os.system(cmd)
	os.remove("tmp.fa")
def homer_findMotifsGenome(name, optional):
	cmd = "findMotifsGenome.pl tmp_homer.bed hg19 %s %s" % (name, optional)
	os.system(cmd)
	os.remove("tmp_homer.bed")
def amd(name, optional, motif_db):
	cmd = "AMD.bin %s -F tmp.fa -B /c/wanghw/software/AMD-motifjournal.pone.0024576.s004/Bgresult1000.txt " % optional
	print cmd
	os.system(cmd)
	nline = ossystemresult("cat tmp.fa|wc -l")
	nline = int(nline)
	os.rename("tmp.fa", "%s.fa" % name)
	os.rename("tmp.fa.Matrix", "%s.Matrix" % name)
	os.rename("tmp.fa.Details", "%s.Details" % name)
	matrixfp = open("%s.Matrix" % name)
	motifmemefp= open("%s.meme" % name, "w")
	motifmatrix={}
	nlinematrix={}
	switch=False
	for line in matrixfp:
		if line.rstrip().endswith(":"):
			a=line.rstrip().rstrip(":")
			nlinematrix[a]=0
			motifmatrix[a]=""
			switch=True
		elif switch == True:
			motifmatrix[a]+="\t".join(line.split("\t")[1:])
			nlinematrix[a]+=1
	print motifmatrix
	print nlinematrix
	modelhead='''MEME version 4.4

ALPHABET= ACGT

strands: + -

Background letter frequencies (from uniform background):
A 0.25000 C 0.25000 G 0.25000 T 0.25000
'''
	modelmotif='''MOTIF %s %s

letter-probability matrix: alength= 4 w= %d nsites= %d E= 0
%s
'''
	print >> motifmemefp, modelhead
	for motifname in motifmatrix.keys():
		print >> motifmemefp, modelmotif % (motifname, motifname, nlinematrix[motifname], nline, motifmatrix[motifname])
	motifmemefp.close()
	matrixfp.close()
	cmd_tomtom = "tomtom -o %s %s %s" % (("%s_tomtom_out" % name),("%s.meme" % name),motif_db)
	os.system(cmd_tomtom)


if __name__ == '__main__':
	starttime=time.time()
	if len(sys.argv) < 2:
		print USAGE % sys.argv[0]
		sys.exit(1)
	
	# defaults
	motif_db = "/c/wanghw/motif_database/Transfect_9.2.meme "
	sequence = "/d/database/hg19/hg19.fa.masked"
	optional_database={"homer": " -mask -S 20 -len 8,10,12,14 -p 4 -size 200", "meme-chip": " -meme-nmotifs 5 -meme-minw 6 -meme-maxw 20", "amd": "-T 50 -CO 0.6 -FC 1.2"}
	tools = "meme-chip"
	d=200
	
	opts,args=getopt.getopt(sys.argv[1:],"",["tools=","motif_db=","sequence=","output=","d="])
	for o,a in opts:
		if o == '--tools':
			tools = a
			optional = optional_database[tools]
			name = os.path.basename(args[0]).replace("_p_s.bed","")+"_"+tools+"_output"
		elif o=='--motif_db':
			motif_db=a
		elif o == '--sequence' :
			sequence = a
		elif o=='--output':
			name  =a +"_"+tools+"_output"
		elif o=='--d' :
			d = a
	psfp = open(args[0])
	if len(args) !=1 :
		optional = " ".join(args[1:])
		
	if tools == "meme-chip":
		d = int(d)
		getfasta_centor(psfp, d, sequence)
		meme_chip(name, motif_db, optional)
	elif tools == "homer":
		d = int(d)
		getbed_centor(psfp, d)
		homer_findMotifsGenome(name, optional)
	elif tools == "amd":
		if d == "fullength":
		    getfasta_fulllength(psfp, sequence)
		else:
		    d = int(d)
		    getfasta_centor(psfp, d, sequence)
		amd(name, optional, motif_db)
		
	endtime=time.time()
	print "it takes %d seconds to run this program!" % (endtime-starttime)
