#!/bin/bash
#sh build-index.sh [path to genome sequence fasta file]  [path to gene annotation gtf file]
#NOTE: the chromsome name of genome sequence file and gene annotation file must be same. It can be "chr1" or "1".
#zjuwhw:zju_whw@gmail.com
#2015-06-08
DNAREF=$1 #absolute path
GTF=$2 #absolute path

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
