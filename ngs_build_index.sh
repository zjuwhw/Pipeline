#!/bin/bash
#sh build-index.sh [path to genome sequence fasta file]  [path to gene annotation gtf file]
#NOTE: the chromsome name of genome sequence file and gene annotation file must be same. It can be "chr1" or "1".
#zjuwhw:zju_whw@gmail.com
#2015-06-08
alias awk="awk -F '\t' -v OFS='\t'"

DNAREF=$1 #absolute path to genome sequence fasta file
GTF=$2 #absolute path to gene annotation gtf file

dir=$(dirname $DNAREF)
basename=$(echo $DNAREF|awk '{split($1,x,"/");print x[length(x)]}')
gtfbasename=$(echo $GTF|awk '{split($1,x,"/");print x[length(x)]}')

#bwa index
echo "building bwa index ..."
mkdir ${dir}/bwaindex
ln -s $DNAREF ${dir}/bwaindex/
bwa index -a bwtsw ${dir}/bwaindex/${basename}

#bowtie2 index
echo "building bowtie2 index ..."
mkdir ${dir}/bowtie2index
ln -s $DNAREF ${dir}/bowtie2index/
bowtie2-build ${dir}/bowtie2index/${basename} ${dir}/bowtie2index/${basename%.fa}

#tophat transcriptome index
echo "building tophat transcriptome index ..."
mkdir ${dir}/TophatTranscriptomeIndex
tophat -p 6 -G $GTF --transcriptome-index=${dir}/TophatTranscriptomeIndex/${gtfbasename%.gtf}.RNA ${dir}/bowtie2index/${basename%.fa}
rm -rf tophat_out

#STAR index
echo "building STAR index ..."
mkdir ${dir}/STARindex
STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir ${dir}/STARindex/ \
--genomeFastaFiles $DNAREF \
--sjdbGTFfile $GTF \
--sjdbOverhang 100 

#HISAT index
echo "building HISAT index ..."
mkdir ${dir}/hisatindex
ln -s $DNAREF ${dir}/hisatindex/
hisat-build ${dir}/hisatindex/${basename} ${dir}/hisatindex/${basename%.fa}

#fasta index for samtools and bedtools
samtools faidx $DNAREF

#chrom size file
cut -f1-2 ${DNAREF}.fai > ${DNAREF}.chrom.sizes

#fasta index for picard
samtools dict $DNAREF > ${DNAREF}.dict
#picard CreateSequenceDictionary R=$DNAREF O=${DNAREF}.dict

#refFlat for picard.CollectRnaSeqMetrics using gtfToGenePred
gtfToGenePred -ignoreGroupsWithoutExons $GTF ${GTF%gtf}genePred
awk '{print $1, $0}' ${GTF%gtf}genePred > ${GTF%gtf}refFlat
genePredToBed ${GTF%gtf}genePred ${GTF%gtf}bed

#gc content of transcript sequence fasta
bedtools getfasta -name -split -fi ${DNAREF} -bed ${GTF%gtf}bed -fo ${DNAREF}.rna.fasta
python calculate_cg_fasta.py ${DNAREF}.rna.fasta


#rRNA interval list for picard.CollectRnaSeqMetrics using code from
#https://gist.github.com/slowkow/b11c28796508f03cdf4b; https://www.biostars.org/p/120145/
perl -lane 'print "\@SQ\tSN:$F[0]\tLN:$F[1]\tAS:hg38"' ${DNAREF}.chrom.sizes > ${GTF%gtf}.rRNA.interval_list

grep 'gene_biotype "rRNA"' $GTF |\
awk '$3=="transcript"' |\
cut -f1,4,5,7,9 |\
perl -lane '/transcript_id "([^"]+)"/ or die "no transcript_id on $."; print join "\t", (@F[0,1,2,3], $1)'|\
sort -k1V -k2n -k3n  >> ${GTF%gtf}.rRNA.interval_list
