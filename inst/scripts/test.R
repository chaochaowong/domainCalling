#'
#' test do_bowtie2_hg38
#'
pkgDir <- "~/tapscott/CutAndRun/hg19.CutandRun2"
ngsDir <- "/fh/fast/henikoff_s/SR/ngs/illumina/henikoff/170622_SN367_0947_BHN3VMBCXY/Unaligned/Project_henikoff"
fqDir <- list.files(ngsDir, pattern="Sample_RR", full.name=TRUE)
sampleDir  <-  fqDir[4]
cmd <- sprintf("sbatch -n6 ~/tapscott/R_package/CutAndRunPipelineTools/inst/scripts/do_bowtie2_hg38.sh %s %s", sampleDir, pkgDir)
system(cmd)



#'
#' get sample information
#'
pkgDir <- "/fh/fast/tapscott_s/CompBio/CutAndRun/hg19.CutAndRun3.XY"
rresnickDir <- "/fh/fast/tapscott_s/CompBio/rresnick/CutandRun"
bamDir <- file.path(rresnickDir, c("0525_dig1",  "0808_dig2"))

bamFiles1 <- file.path(bamDir[1], c("resequence/XYrenorm/RR_HsDm_XY_0525.bam",
                                   "resequence/XYrenorm/RR_HsDm_XY_0525.sam.dm6",
                                   "Analysis_0525_Dig1/RR_HsDm_2o_0525.bam",
                                    "Analysis_0525_Dig1/RR_HsDm_2o_0525.sam.dm6"))

bamFiles2 <- file.path(bamDir[2], c("Analysis_0808_Dig2/RR_HsDm_XY1.bam",
                                     "Analysis_0808_Dig2/RR_HsDm_XY1.sam.dm6",
                       "Analysis_0808_Dig2/RR_HsDm_XY2.bam",
                       "Analysis_0808_Dig2/RR_HsDm_XY2.sam.dm6",
                       "Analysis_0808_Dig2/RR_HsDm_2.bam",
                       "Analysis_0808_Dig2/RR_HsDm_2.sam.dm6"))

si <- data.frame(sample_name=c(basename(bamFiles1), basename(bamFiles2)),
                 file_bam=c(bamFiles1, bamFiles2),
                 genome=rep(c("hg19", "dm6"), 5),
                 pheno_type=c("XY", "XY", "neg", "neg",
                     "XY", "XY", "XY", "XY", "neg", "neg"),
                 prep_date=c(rep("052517", 4), rep("080817", 6)),
                 stringsAsFactors=FALSE)

#'
#' test makeBedGraphPerSample()
#' 
makeBedGraphPerSample(bam_file, scaling_factor=NULL,
                      bigWig=TRUE,
                      remove.duplicate=TRUE,
                      outDir=file.path(pkgDir, "coverage"))

makeBedGraphPerSample(bam_file,
                      scaling_factor=0.9,
                      bigWig=FALSE,
                      remove.duplicate=FALSE,
                      outDir=file.path(pkgDir, "coverage"))
#'
#' test makeBedGraph()
#' 
makeBedGraphPerSample(bam_file, scaling_factor=NULL,
                      bigWig=FALSE,
                      remove.duplicate=FALSE,
                      outDir=file.path(pkgDir, "coverage"))

#'
#' (1) simple peak calling
#'

#'
#' use csaw to mimic peak calling
#' 
