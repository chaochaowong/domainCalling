#!/bin/bash
#./do_bowtie2_hg38.sh
#SBATCH -n6 -t 1-0 -p campus --mail-type=ALL --mail-user=cwon2@fhcrc.org -A tapscott_s

sampleDir=$1
sampleName=`basename ${sampleDir}`
pkgDir=$2
bowtie2OutDir=$pkgDir/bowtie2_hg38
fastqcDir=$pkgDir/fastqc

echo $sampleName
mkdir -p $bowtie2OutDir
mkdir -p $fastqcDir

hg38_genomeBuild=/shared/biodata/ngs/Reference/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome
dm6_genomeBuild=/shared/biodata/ngs/Reference/iGenomes/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome
bowtieVersion=/home/solexa/apps/bowtie/bowtie2-2.2.6
tophatVersion=/home/solexa/apps/tophat/tophat-2.1.0.Linux_x86_64

#umask 0002
export PATH=/usr/kerberos/bin:/usr/local/bin:/bin:/usr/bin:/usr/X11R6/bin:/opt/moab/bin:$bowtieVersion:$tophatVersion:/home/solexa/apps/samtools/samtools-0.1.19:/home/solexa/apps/FastQC
#remove /app/bin from the PATH, as it will default to running bowtie 0.12.8
export PATH=${PATH/\/app\/bin:}

# zcat all the files
cd $sampleDir

fq1=$bowtie2OutDir/$sampleName\_R1.fastq
fq2=$bowtie2OutDir/$sampleName\_R2.fastq

zcat $(ls *_R1_*.fastq.gz) > $fq1
zcat $(ls *_R2_*.fastq.gz) > $fq2

#align to hg38
cd $bowtie2OutDir

bowtie2 --local --very-sensitive-local --no-unal --no-mixed --no-discordant -q --phred33 -I 10 -X 700 --threads 6 -x $hg38_genomeBuild -1 $fq1 -2 $fq2 > $sampleName.sam

# sort and index
samtools view -Sb $sampleName.sam > $sampleName.bam
rm $sampleName.sam

samtools sort -@ 4 $sampleName.bam $sampleName.bam.sorted
mv $sampleName.bam.sorted.bam $sampleName.bam
samtools index $sampleName.bam
touch $sampleName.bowtie2Done.txt

#align to dm6
spikeName=$sampleName.dm6
bowtie2 --local --very-sensitive-local --no-unal --no-mixed --no-discordant -q --phred33 -I 10 -X 700 --threads 6 -x $dm6_genomeBuild -1 $fq1 -2 $fq2 > $spikeName.sam

samtools view -Sb $spikeName.sam > $spikeName.bam
rm $spikeName.sam
rm $fq1
rm $fq2

# sort and index dm6 bam
samtools sort -@ 4 $spikeName.bam $spikeName.bam.sorted
mv $spikeName.bam.sorted.bam $spikeName.bam
samtools index $spikeName.bam
touch $spikeName.bowtie2Done.txt

exit 0

