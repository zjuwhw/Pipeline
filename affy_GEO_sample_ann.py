USAGE = '''GEO_sample_ann.py --- extract the sample information from GEO_series_matrix.txt file
USAGE:
    python %s GEO_series_matrix.txt
    
'''
import os, sys, time, re

def allsame(a):
    tmp = a[0]
    switch = True
    for i in range(len(a)):
        if a[i] != tmp:
            switch = False
    return switch
def unquote(string):
    return string.strip('"')

if __name__ == '__main__':
    starttime=time.time()
    if len(sys.argv) < 2:
        print USAGE % sys.argv[0]
        sys.exit(1)
    
    f = open(sys.argv[1])
    f2 = open(sys.argv[1].replace("series_matrix.txt","sample_matrix.txt"),"w")
    f3 = open(sys.argv[1].replace("series_matrix.txt","value_matrix.txt"),"w")
    container1 = []
    for line in f:
        if line.startswith("!Sample_"):
            linelist = line.lstrip("!").rstrip().split("\t")
            if not allsame(linelist[1:]) and (not linelist[0].startswith("Sample_contact")) and (not linelist[0].endswith("file")):
                linelist = map(unquote, linelist)
                container1.append(linelist)
        elif not line.startswith("!"):
            linelist = line.rstrip().split("\t")
            linelist = map(unquote, linelist)
            if linelist != [""]:
                print >> f3, "\t".join(linelist)
    tcontainer1 = map(list, zip(*container1))
    for newline in tcontainer1:
        print >> f2, "\t".join(newline)
    f.close()
    f2.close()
    f3.close()
    
    


