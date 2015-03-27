#author: zjuwhw
#date: 2015-03-27
usage='''chipseq_trans2gene.py --- input _p_s.bed file and output the peak2gene file

USAGE:
    python %s [--type=#] [--promoter_up=#] [--promoter_dn=#] [--intergenic_range=#] [--GeneAnnBed12=<file>] [--output=<file>] psbedfile

#the input file is the _p_s.bed peak bed file
#--type: there are three types to associate the peak with gene, "peak2gene", "gene2peak", "peakAroundgene". Default: peak2gene
#--promoter_up: the upstream range of promoter defination. Default: -2000
#--promoter_dn: the dnstream range of promoter defination. Default: 2000
#--intergenic_range: the both range of the intergenic region. Default: 20000
#--GeneAnnBed12: the gene annotation bed12 file. Default: /c/wanghw/annotation/refseq_hg19_07292013.bed
#--output: the output file name. Default: the basename of psbedfile + type

Note:
peak2gene type: every peak (including the intergenic peak) will be retained once, but gene maybe exist multiple times.
gene2peak type: every gene will be retianed once, and one peak maybe exist multiple times, some peaks will be removed, peak can be in intergenic region.
peakAroundgene type: Both gene and peak can exist multiple times, and all intergenic peak will be removed.

Default:
    python %s --type peak2gene --promoter_up -2000 --promoter_dn 2000 --intergenic_range 20000 --GeneAnnBed12 /c/wanghw/annotation/refseq_hg19_07292013.bed _p_s.bed
'''

import os,sys,getopt,time

def ps2summit(psfile):
    f=open(psfile)
    f2=open(os.path.basename(psfile).replace("_p_s.bed","_summit.bed6"),"w")
    for line in f:
        linelist=line.rstrip().split("\t")
        print >>f2, "%s\t%s\t%s\t%s\t%s\t%s" %(linelist[0],str(int(linelist[4])-1),linelist[4],linelist[3],linelist[1],linelist[2])
    f.close()
    f2.close()
    return os.path.basename(psfile).replace("_p_s.bed","_summit.bed6")
def GeneAnnBed12ToTSSbed14(geneannfile):
    f=open(geneannfile)
    f2=open(os.path.basename(geneannfile).replace(".bed",".TSS.bed14"),"w")
    for line in f:
        linelist = line.rstrip().split("\t")
        linelist.append(linelist[1])
        linelist.append(linelist[2])
        if linelist[5] == "+":            
            linelist[2] = str(int(linelist[1]) + 1)
        else:
            linelist[1] = str(int(linelist[2]) - 1)
        print >> f2, "\t".join(linelist)
    f.close()
    f2.close()
    return os.path.basename(geneannfile).replace(".bed",".TSS.bed14")
def peak2gene(peakfile_summit_bed6, genefile_TSS_bed14, peak2gene_name):
    cmd = r'''bedtools closest -D b -t first -a %s -b %s |awk -F '\t' -v OFS='\t' '$1!="." && $7!="."{print $0}' |sort -k1,1 -k2,2n -k3,3n -k 8,8n -k9,9n > %s ''' % (peakfile_summit_bed6, genefile_TSS_bed14, peak2gene_name)
    os.system(cmd)
    return peak2gene_name
def gene2peak(peakfile_summit_bed6, genefile_TSS_bed14, gene2peak_name):
    cmd = r'''bedtools closest -D a -t first -a %s -b %s |awk -F '\t' -v OFS='\t' '{print $15,$16,$17,$18,$19,$20,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$21}' |awk -F '\t' -v OFS='\t' '$1!="." && $7!="."{print $0}' |sort -k1,1 -k2,2n -k3,3n -k 8,8n -k9,9n  > %s ''' % (genefile_TSS_bed14, peakfile_summit_bed6, gene2peak_name)
    os.system(cmd)
    return gene2peak_name
def peakAroundgene(peakfile_summit_bed6, genefile_TSS_bed14, window, peakAroundgene_name):
    cmd = r'''bedtools window -w %s -a %s -b %s| awk -F '\t' -v OFS='\t' '$1!="." && $7!="."{a=$3-$9;if($12=="-"){a=-a}else{a=a};print $0,a}' |sort -k1,1 -k2,2n -k3,3n -k 8,8n -k9,9n > %s ''' % (window, peakfile_summit_bed6, genefile_TSS_bed14, peakAroundgene_name)
    os.system(cmd)
    return peakAroundgene_name
def summit_in_exon(d, string_exon_length, string_exon_str):
    result = False
    exon_length_list = string_exon_length.rstrip(",").split(",")
    exon_str_list = string_exon_str.rstrip(",").split(",")
    for i in range(len(exon_length_list)):
        if d > int(exon_str_list[i]) and d <= int(exon_str_list[i])+int(exon_length_list[i]):
            result = True
    return result
def summit_in_intron(d, string_exon_length, string_exon_str):
    result = False
    exon_length_list = string_exon_length.rstrip(",").split(",")
    exon_str_list = string_exon_str.rstrip(",").split(",")
    if len(exon_length_list) != 1:
        for i in range(len(exon_length_list)-1):
            if d > int(exon_str_list[i]) + int(exon_length_list[i]) and d <= int(exon_str_list[i+1]):
                result = True
    return result
    
def add_information(tmp_bed21, outfilename):
    f = open(tmp_bed21)
    f2 = open(outfilename, "w")
    print >> f2, "\t".join(["#peak_chr","peak_str","peak_end","peak_id","peak_summit","gene_chr","gene_str","gene_end","gene_symbol","gene_refseqid","gene_strand","summit2tss","position_type","multipeak_count","multipeak_index","multigene_count","multigene_index","multitrans_count","multitrans_index"])
    multipeak_count = []
    multipeak_index1 = {}
    multipeak_index2 = {}
    multigene_count = []
    multigene_index1 = {}
    multigene_index2 = {}
    multitrans_count = []
    multitrans_index1 = {}
    multitrans_index2 = {}
    for line in f:
        linelist = line.rstrip().split("\t")
        multipeak = linelist[3]
        multipeak_index1[multipeak] = 0
        multipeak_index2[multipeak] = 0
        multigene = linelist[9]
        multigene_index1[multigene] = 0
        multigene_index2[multigene] = 0
        multitrans = linelist[10]
        multitrans_index1[multitrans] = 0
        multitrans_index2[multitrans] = 0
    f.seek(0)
    for line in f:
        linelist = line.rstrip().split("\t")
        multipeak = linelist[3]
        if multipeak not in multipeak_count:
            multipeak_count.append(multipeak)
        multipeak_index2[multipeak] += 1
        multigene = linelist[9]
        if multigene not in multigene_count:
            multigene_count.append(multigene)
        multigene_index2[multigene] += 1
        multitrans = linelist[10]
        if multitrans not in multitrans_count:
            multitrans_count.append(multitrans)
        multitrans_index2[multitrans] += 1
    f.seek(0)
    for line in f:
        linelist = line.rstrip().split("\t")
        #position_type
        if int(linelist[20]) >= promoter_up and int(linelist[20]) <= promoter_dn:
            position_type = "promoter"
        elif summit_in_exon(int(linelist[2])-int(linelist[18]), linelist[16], linelist[17]):
            position_type = "exon"
        elif summit_in_intron(int(linelist[2])-int(linelist[18]), linelist[16], linelist[17]):
            position_type = "intron"
        elif int(linelist[20]) < promoter_up and int(linelist[20]) >= -intergenic_range:
            position_type = "upstream"
        elif int(linelist[2]) > int(linelist[18])-intergenic_range and int(linelist[2]) <= int(linelist[19]) + intergenic_range:
            position_type = "downstream"
        else:
            position_type = "intergenic"
        #multipeak_count, multipeak_index, multigene_count, multigene_index, multitrans_count, multitrans_index
        multipeak = linelist[3]
        multipeak_index1[multipeak] += 1
        multigene = linelist[9]
        multigene_index1[multigene] += 1
        multitrans = linelist[10]
        multitrans_index1[multitrans] += 1
        outputline = [linelist[0], linelist[4], linelist[5], linelist[3], linelist[2], linelist[6], linelist[12], linelist[13], linelist[9], linelist[10], linelist[11], linelist[20], position_type, "%s" % (multipeak_count.index(multipeak)+1), "%sof%s" % (multipeak_index1[multipeak], multigene_index2[multigene]), "%s" % (multigene_count.index(multigene)+1), "%sof%s" % (multigene_index1[multigene], multigene_index2[multigene]), "%s" % (multitrans_count.index(multitrans)+1), "%sof%s" % (multitrans_index1[multitrans], multitrans_index2[multitrans])]
        print >> f2, "\t".join(outputline)
    f.close()
    f2.close()
    
if __name__ == '__main__':
    
    starttime=time.time()
    if len(sys.argv) < 2:
        print usage % (sys.argv[0],sys.argv[0])
        sys.exit(1)
        
    opts, args = getopt.getopt(sys.argv[1:], "", ["type=","promoter_up=","promoter_dn=","intergenic_range=","GeneAnnBed12=","output="] )
    inputfile=args[0]
    
    # defaults
    types=["peak2gene","gene2peak","peakAroundgene"]
    ntype="peak2gene"
    promoter_up=-2000
    promoter_dn=2000
    intergenic_range=20000
    GeneAnn="/c/wanghw/annotation/refseq_hg19_07292013.bed"
    outfilename=os.path.basename(inputfile).replace("_p_s.bed","") + "_" + ntype + ".txt"
    
    for o,a in opts:
        if o == '--type':
            if a in types:
                ntype = a
                outfilename=os.path.basename(inputfile).replace("_p_s.bed","") + "_" + ntype + ".txt"
            else:
                print 'The type must be one of "peak2gene", "gene2peak", "peakAroundgene"'
                sys.exit(1)
        elif o == '--promoter_up':
            promoter_up = int(a)
        elif o == '--promoter_dn':
            promoter_dn = int(a)
        elif o == '--intergenic_range':
            intergenic_range = int(a)
        elif o == '--GeneAnnBed12':
            GeneAnn = a
        elif o == '--output':
            outfilename = a
    
    intermediates=[]
    intermediates.append(ps2summit(inputfile))
    intermediates.append(GeneAnnBed12ToTSSbed14(GeneAnn))
    if ntype == "peak2gene":
        intermediates.append(peak2gene(intermediates[0], intermediates[1],os.path.basename(inputfile).replace("_p_s.bed","")+"_tmp.peak2gene.bed21"))
    elif ntype == "gene2peak":
        intermediates.append(gene2peak(intermediates[0], intermediates[1],os.path.basename(inputfile).replace("_p_s.bed","")+"_tmp.gene2peak.bed21"))
    elif ntype == "peakAroundgene":
        intermediates.append(peakAroundgene(intermediates[0], intermediates[1], intergenic_range, os.path.basename(inputfile).replace("_p_s.bed","")+"_tmp.peakAroundgene.bed21"))
    add_information(intermediates[2],outfilename)
    for intermediate in intermediates:
        os.remove(intermediate)
    endtime = time.time()
    print "it takes %.3f secondes or %.3f minutes" % (endtime -starttime, (endtime-starttime)/60)
    
