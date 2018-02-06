#!/bin/bash
#./do_bowtie2_FLY.sh

fq1=$1
fq2=$2
bowtie2OutDir=$3

base=`basename ${fq1}`
sampleName=`echo ${base//_R1.fastq}`
echo $sampleName

gtfFile=/shared/biodata/ngs/Reference/iGenomes/Drosophila_melanogaster/UCSC/dm6/Annotation/Archives/archive-2015-07-24-09-25-49/Genes/genes.gtf
genomeBuild=/shared/biodata/ngs/Reference/iGenomes/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome
bowtieVersion=/home/solexa/apps/bowtie/bowtie2-2.2.6
tophatVersion=/home/solexa/apps/tophat/tophat-2.1.0.Linux_x86_64

#umask 0002
export PATH=/usr/kerberos/bin:/usr/local/bin:/bin:/usr/bin:/usr/X11R6/bin:/opt/moab/bin:$bowtieVersion:$tophatVersion:/home/solexa/apps/samtools/samtools-0.1.19:/home/solexa/apps/FastQC
#remove /app/bin from the PATH, as it will default to running bowtie 0.12.8
export PATH=${PATH/\/app\/bin:}

#alignment using bowtie2
cd $bowtie2OutDir

bowtie2 --local --very-sensitive-local --no-unal --no-mixed --no-discordant -q --phred33 -I 10 -X 700 --threads 12 -x $genomeBuild -1 $fq1 -2 $fq2 > $sampleName.sam

samtools view -Sb $sampleName.sam > $sampleName.bam.FLY
rm $sampleName.sam

#samtools sort -@ 4 $sampleName.bam $sampleName.bam.sorted
#mv $sampleName.bam.sorted.bam $sampleName.bam
#samtools index $sampleName.bam
touch $sampleName.bowtie2Done.txt
exit 0

