# domainCalling

_domainCalling_ is a R/Bioconductor package, and it has been tested on R 3.4.0/Bioc 3.5. This package is designed for
detection of broad range of non-coding transcripts (from total RNA-seq) or binding regions, particularly from CutAndRun-seq (CAR-seq) data.

### Requirement
You need the _csaw_, _edgeR_ and _ChIPseeker_ packages. A TxDb and a compatible org.xx.db package for annotation.

# Using _csawDomainCalling()_
This function uses window-based counting method provided by the _csaw_ package. It first counts the reads overlapping with a 
sliding window. The second step is to filter out uninteresting windows in which the abounance of the non-background samples is
lower then 3 fold change to the backgroud or does not exceed the threshold. Finally, the retained windows are merged with 
neighbors within certain perimeters. The merged regions is the potential bind regions of the protein.

## Sample Information _domainCalling::getSampleInfo()_

## Spike factor _domainCalling:::addSpikeFactor()_

## Workflow
