#'
#' Example of using do_bowtie2_hg38.sh to do alignments
#' 
#' Requirement: (1) The sample directory that store FASTQ files of a simple.
#'              (2) The main direcotory of your analaysis.
#' Return products: Create bowtie2_hgxx sub-directory and yield sampleName.bam and
#'                   sampleName.dm6.bam (reads align to fly).
#' Recommendation: use sbatch to run the script.
#' 

pkgDir <- "~/tapscott/CutAndRun/hg19.CutandRun2"
ngsDir <- "/fh/fast/henikoff_s/SR/ngs/illumina/henikoff/170622_SN367_0947_BHN3VMBCXY/Unaligned/Project_henikoff"
fqDir <- list.files(ngsDir, pattern="Sample_RR", full.name=TRUE)
sampleDir  <-  fqDir[4]
cmd <- sprintf("sbatch -n6 ~/tapscott/R_package/CutAndRunPipelineTools/inst/scripts/do_bowtie2_hg38.sh %s %s", sampleDir, pkgDir)
system(cmd)
