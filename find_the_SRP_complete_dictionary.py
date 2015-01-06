USAGE='''find_the_SRP_complete_dictionary.py --- input a SRP number, output the complete dictionary of all sra files in the SRP number folder

USAGE:
    python %s SRPnumber(SRPXXXXXX)

NOTES:
the SRPnumber format must be like "SRPXXXXXX"(SRP + 6 numbers)
'''
import os,sys,getopt,time,re,ftplib

def ossystemresult(command):
    fp=os.popen(command,"r")
    return fp.read()

		
def listdictionary(link):
	nlink = ossystemresult("curl -l "  + link.strip("/") + "/" )
	return nlink.strip().split("\n")

def isdictionary(link):
	nlinklist = listdictionary(link)
	if len(nlinklist) == 1 and nlinklist[0] == link.strip("/").split("/")[-1]:
		return False
	else:
		return True

def search_dic(link) :
	link = link.strip("/") + "/"
	m = listdictionary(link)
	for i in m:
		print >> f, link + i
		if not i.endswith("sra"):
			search_dic(link + i)
			
if __name__ == '__main__':
	starttime=time.time()
	if len(sys.argv) < 2:
		print USAGE % sys.argv[0]
		sys.exit(1)
    #opts, args = getopt.getopt(sys.argv[1:], "", ["type=","strand=","nbins=","range=","output="])        
    
	match = re.match(r'^SRP[0-9]{6}$', sys.argv[1])
	if not match :
		print '''the SRPnumber %s is not like "SRPXXXXXX"(SRP + 6 numbers) format''' % sys.argv[1]
	SRPnumber = sys.argv[1]
	SRPnumber3 = sys.argv[1][3:6]
	SRPlink = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP" + SRPnumber3 + "/" + SRPnumber
	global f
	f = open(sys.argv[1]+"_list.txt","w")
	search_dic(SRPlink)
	f.close()
	endtime = time.time()
	print starttime, endtime, endtime -starttime
