#' This is an example of running the pipeline.
#' Fist, assuming sorted and indexed BAM files of the samples and spike sample are
#' provided (bowtie2_hgxx). 
#' ml R/3.4.0-foss-2016b-fh1
library(parallel)
library(Rsamtools)
source("/fh/fast/tapscott_s/CompBio/R_package/CutAndRunPipelineTools/R/makeBedGraphPerSample.R")
source("/fh/fast/tapscott_s/CompBio/R_package/CutAndRunPipelineTools/R/getSampleInformation.R")
source("/fh/fast/tapscott_s/CompBio/R_package/CutAndRunPipelineTools/R/peakCallingFromBed.R")
source("/fh/fast/tapscott_s/CompBio/R_package/CutAndRunPipelineTools/R/tools.R")
#'
#' (1) Sample Information
#'     Create a data.frame with basic information: sample_name, file_bam and spike_bam.
#'     Call up getSampleInformation
#' 
pkgDir <- "/fh/fast/tapscott_s/CompBio/CutAndRun/hg19.CutandRun2"
bamFiles <- list.files(file.path(pkgDir, "bowtie2_hg19"), pattern="\\.bam$",
                       full.name=TRUE)
names(bamFiles) <- sub(".bam", "", basename(bamFiles))
spikeFiles <- list.files(file.path(pkgDir, "bowtie2_dm6"), pattern="\\.bam.FLY$",
                         full.name=TRUE)
names(spikeFiles) <- sub(".bam.FLY", "", basename(spikeFiles))
si <- data.frame(sample_name=names(bamFiles),
                 file_bam=bamFiles,
                 spike_bam=spikeFiles[names(bamFiles)],
                 group=c("neg", "Dux4", "K27Ac", "XY"),
                 stringsAsFactors=FALSE)


si <- getSampleInfo(sample_info=si, multipling_factor=21000, cores=4L)
sample_info <- si
save(sample_info, file=file.path(pkgDir, "data", "sample_info.rda"))

#'
#' Use edgeR to determine the norm_factor
#'
norm_factors <- getSpikeNormFactor(sample_info, fragment=180)

#'
#' (2) make BedGraph: one by one or use the wrapper function makeBedGraph
#' 
source("/fh/fast/tapscott_s/CompBio/R_package/CutAndRunPipelineTools/R/makeBedGraphPerSample.R")
bdgDir <- file.path(pkgDir, "bedGraph")

#' makeBedGraphPerSample(): do it one-by-one
makeBedGraphPerSample(si$file_bam[4], scaling_factor=NULL, bigWig=TRUE,
                      remove.duplicate=FALSE, outDir=bdgDir)
makeBedGraphPerSample(si$file_bam[4], scaling_factor=0.9, bigWig=TRUE,
                      remove.duplicate=TRUE, outDir=bdgDir)

#' makeBedGraph(): if it is not working that means each core does not enough memory to
#' support. Use lapply and makeBedGraphPerSample instead or
#' use grabnode, select more then four cores, and allocate 8 GB to each core
makeBedGraph(sample_info=si, scale=FALSE, remove.duplicate=FALSE, outDir=bdgDir,
             bigWig=TRUE, cores=4L) # not working when I allocate 8GB to each core

#' if "makeBedGraph()" not working due to memory allocation, use lapply or mapply
lapply(sample_info$file_bam, makeBedGraphPerSample,
       scaling_factor=NULL,
       bigWig=FALSE,
       remove.duplicate=FALSE,
       outDir=bdgDir)
#' scalling; use mapply and makeBedGraphPerSample()
mapply(makeBedGraphPerSample,sample_info$file_bam, 
       scaling_factor=sample_info$spike_factor,
       MoreArgs=list(bigWig=FALSE, remove.duplicate=FALSE, outDir=bdgDir))

#'
#' (3) simple peak calling using threshold
#'
source("/fh/fast/tapscott_s/CompBio/R_package/CutAndRunPipelineTools/R/peakCallingFromBed.R")

bdgDir <- file.path(pkgDir, "bedGraph")
bdg_files <- list.files(bdgDir, pattern=".bdg$",
                        full.name=TRUE)
peakDir <- file.path(pkgDir, "peaks") 
bdg_file <- bdg_files[4]
peaks.gr <- peakCallingFromBed(bed_file=bdg_file, format="bedGraph",
                               threshold=10L, standard.chromosome=TRUE,
                               species="Homo_sapiens",
                               min.range=100L, max.range=20000,
                               gapwidth=2000L,
                               outDir=peakDir)


#'
#' (4) Annotation
#'
library(ChIPseeker)
source("/fh/fast/tapscott_s/CompBio/R_package/CutAndRunPipelineTools/R/peakAnnotationEnbs.R")
library(EnsDb.Hsapiens.v75)
seqlevelsStyle(EnsDb.Hsapiens.v75) <- "UCSC"
library(TxDb.Hsapiens.UCSC.hg19.ensGene)
seqlevelsStyle(TxDb.Hsapiens.UCSC.hg19.ensGene) <- "UCSC"
load("~/tapscott/hg19/hg19.rmsk.rda")

#' NOTE: At this moment, both TxDb and EnsDb are both needed. Once
#' ChIPseeker releases newer version, then TxDb will no longer being needed.
prefix <- "Sample_RR_HsDm_XY_0525_scalled0.9_peaks_thres10"
peaks <- peakAnnotationFromGR(peaks.gr, TxDb=TxDb.Hsapiens.UCSC.hg19.ensGene,
                              EnsDb=EnsDb.Hsapiens.v75,
                              rmsk=hg19.rmsk, annotate.repeat=TRUE)
plotAnnoBar(peaks$peakAnno)
vennpie(peaks$peakAnno)
write.csv(as.data.frame(peaks$peaks.gr),
          file=file.path(peakDir, "Sample_RR_HsDm_XY_0525_scalled0.9_peaks_thres10.csv"))
save(peaks, file=file.path(peakDir, "Sample_RR_HsDm_XY_0525_scalled0.9_peaks_thres10.rda"))

#'
#' (5) peak calling using csaw pipeline
#' 
source("/fh/fast/tapscott_s/CompBio/R_package/CutAndRunPipelineTools/R/peakCalling.R")
source("/fh/fast/tapscott_s/CompBio/R_package/CutAndRunPipelineTools/R/getSampleInformation.R")
load(file.path(pkgDir, "data", "sample_info.rda"))
cores <- 2L

library(EnsDb.Hsapiens.v75)
seqlevelsStyle(EnsDb.Hsapiens.v75) <- "UCSC"
standard_chromosome <- standardChromosomes(EnsDb.Hsapiens.v75, species="Homo_sapiens")
pkgDir <- "/fh/fast/tapscott_s/CompBio/CutAndRun/hg19.CutandRun2"

#' one sample without spike normalization
si <- sample_info[4, ]
norm_factors <- getSpikeNormFactor(sample_info[c(1,4), ], fragment=180)

domains1 <- domainCallingFromBam(sample_info=si, 
                            background_idx=NULL,
                            threshold=10L,
                            spacing=100, 
                            dedup=FALSE,
                            standard.chromosome=NULL,
                            spikeNorm=FALSE,
                            filter.By.GlobalEnrichment=TRUE,
                            max.width=30000,
                            window.gapwidth=2000)
#' one sample with spike normalization
domains2 <- domainCallingFromBam(sample_info=si, 
                                 background_idx=NULL,
                                 threshold=10L,
                                 spacing=100, 
                                 dedup=FALSE,
                                 standard.chromosome=NULL,
                                 spikeNorm=TRUE,
                                 norm.factors=1/0.9,
                                 filter.By.GlobalEnrichment=TRUE,
                                 max.width=30000,
                                 window.gapwidth=2000)

#' filter based on cluster width > 100? those could just be noise??


#' tesing: no norm.factor. with background sample
si <- sample_info[c(1,4), ]
domains3 <-  csawDomainCalling(sample_info=si, 
                                 background_idx=1,
                                 threshold=10L,
                                 spacing=100, 
                                 dedup=FALSE,
                                 standard.chromosome=NULL,
                                 spikeNorm=FALSE,
                                 filter.By.GlobalEnrichment=FALSE,
                                 max.width=30000,
                                 window.gapwidth=2500)
#' testing using norm.factors and global enrichment
domains4 <-  csawDomainCalling(sample_info=si, 
                                 background_idx=1,
                                 threshold=10L,
                                 spacing=100, 
                                 dedup=FALSE,
                                 standard.chromosome=NULL,
                                 spikeNorm=TRUE,
                                 norm.factors=1/si$spike_factor,
                                 filter.By.GlobalEnrichment=TRUE,
                                 max.width=30000,
                                 window.gapwidth=2500)
