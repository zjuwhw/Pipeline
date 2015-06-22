#inputfasta
import os, sys
inputfasta = sys.argv[1]
output = sys.argv[1].replace(".fa","_gc.txt").replace(".fasta","_gc.txt")
f = open(inputfasta)
f1 = open(output, "w")
switch = True
for line in f:
    if line.startswith(">") and switch:
        switch = False
        header = line.rstrip()[1:]
        string = ""
    elif line.startswith(">") and not switch:
        print >> f1, "%s\t%.8f" % (header, float(string.count("C")+string.count("c")+string.count("G") + string.count("g"))/len(string)  )
        header = line.rstrip()[1:]
        string = ""
    else:
        string += line.rstrip()
print >> f1, "%s\t%.8f" % (header, float(string.count("C")+string.count("c")+string.count("G") + string.count("g"))/len(string)  )
f.close()
f1.close()

