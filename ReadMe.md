# domainCalling

_domainCalling_ is a R/Bioconductor package, and it has been tested on R 3.4.0/Bioc 3.5. This package is designed for
detection of broad range of non-coding transcripts (from total RNA-seq) or binding regions, particularly from CutAndRun-seq (CAR-seq) data.

### Requirement
You need the _csaw_, _edgeR_ and _ChIPseeker_ packages. A TxDb and a compatible org.xx.db package for annotation.

# Using _csawDomainCalling()_
This function uses window-based counting method provided by the _csaw_ package. It first counts reads overlapping with a 
sliding window from non-background and background (if there is any) samples. Secondly it filters out uninteresting windows if
the average abounance of the non-background samples is (1) lower then 3 fold change to the backgroud or (2) does not exceed the 
threshold. Finally, the retained windows are merged with neighbors within designated perimeters. The merged regions is the potential 
binding regions of the protein. 

## Sample Information domainCalling::getSampleInfo()
Providing a data.frame (or DataFrame) containing three columns named `sample_name`, `file_bam`, and `spike_bam` the user can 
use getSampleInfo() to find information from the bam files including fragment size, library size, spike-in counts,
single-ended or pair-ended. These inforamtion is needed to run csawDomainCalling().

Need to add an example here

## Spike factor domainCalling:::getSpikeNormFactor()
The package provides a tool to estimate the normalization factor using spike-in sequence data. 

Need to add an example here

## Workflow
Need to add a workflow here
