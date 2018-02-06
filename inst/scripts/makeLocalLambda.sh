macs2Dir=/fh/fast/tapscott_s/CompBio/CutAndRun/hg19.CutandRun2/MACS2_bam
outDir=/fh/fast/tapscott_s/CompBio/CutAndRun/hg19.CutandRun2/MACS2_bg
testDir=/fh/fast/tapscott_s/CompBio/CutAndRun/hg19.CutandRun2/MACS2_test
bamDir=/fh/fast/tapscott_s/CompBio/CutAndRun/hg19.CutandRun2/bowtie2_hg19

# callpeak from bam files
cd $bamDir
control=Sample_RR_HsDm_2o_0525

sampleName=Sample_RR_HsDm_DUX4_0525
macs2 callpeak -t $sampleName.bam -f BAMPE --outdir $macs2Dir -g hs -n $sampleName -B -q 0.01

sampleName=Sample_RR_HsDm_K27Ac_0525
macs2 callpeak -t $sampleName.bam -f BAMPE --outdir $macs2Dir -g hs -n $sampleName -B -q 0.01 --broad --broad0cutoff 0.05

sampleName=Sample_RR_HsDm_XY_0525
macs2 callpeak -t $sampleName.bam -f BAMPE --outdir $macs2Dir -g hs -n $sampleName -B -q 0.01 --broad --broad-cutoff 0.05

macs2 callpeak -t $sampleName.bam -c $control.bam -f BAMPE --outdir $macs2Dir -g hs -n Sample_RR_HsDm_XY_0525_vsControl -B -q 0.01 --broad --broad-cutoff 0.01


# no scalled: testing
cd $macs2Dir
macs2 bdgcmp -t Sample_RR_HsDm_XY_0525_treat_pileup.bdg -c Sample_RR_HsDm_XY_0525_control_lambda.bdg -m qpois -o $testDir/Sample_RR_HsDm_XY_0525_qvalue.bdg

macs2 bdgbroadcall -i $testDir/Sample_RR_HsDm_XY_0525_qvalue.bdg -c 1.302 -C 1 -l 189 -g 25 -o $testDir/Sample_RR_HsDm_XY_0525_peaks.bed

macs2 bdgcmp -t Sample_RR_HsDm_K27Ac_0525_treat_pileup.bdg -c Sample_RR_HsDm_K27Ac_0525_control_lambda.bdg -m qpois -o $testDir/Sample_RR_HsDm_K27Ac_0525_qvalue.bdg

macs2 bdgbroadcall -i $testDir/Sample_RR_HsDm_K27Ac_0525_qvalue.bdg -c 1.302 -C 1 -l 189 -g 25 -o $testDir/Sample_RR_HsDm_K27Ac_0525_peaks.bed

macs2 bdgcmp -t $outDir/Sample_RR_HsDm_K27Ac_0525_scalled_pileup.bdg -c Sample_RR_HsDm_K27Ac_0525_control_lambda.bdg -m qpois -o $testDir/Sample_RR_HsDm_K27Ac_0525_qvalue2.bdg

macs2 bdgbroadcall -i $testDir/Sample_RR_HsDm_K27Ac_0525_qvalue2.bdg -c 1.302 -C 1 -l 189 -g 25 -o $testDir/Sample_RR_HsDm_K27Ac_0525_peaks2.bed

macs2 bdgcmp -t $outDir/Sample_RR_HsDm_K27Ac_0525_scalled_pileup.bdg -c $outDir/Sample_RR_HsDm_K27Ac_0525_scalled_control_lambda.bdg -m qpois -o $testDir/Sample_RR_HsDm_K27Ac_0525_qvalue3.bdg

macs2 bdgbroadcall -i $testDir/Sample_RR_HsDm_K27Ac_0525_qvalue3.bdg -c 1.302 -C 1 -l 189 -g 25 -o $testDir/Sample_RR_HsDm_K27Ac_0525_peaks3.bed

# call peaks for scalled pileups

#XY
cd $outDir
macs2 bdgcmp -t Sample_RR_HsDm_XY_0525_scalled_pileup.bdg -c Sample_RR_HsDm_XY_0525_scalled_control_lambda.bdg -m qpois -o Sample_RR_HsDm_XY_0525_qvalue.bdg

macs2 bdgbroadcall -i Sample_RR_HsDm_XY_0525_qvalue.bdg -c 1.302 -C 1 -l 182 -g 25 -o Sample_RR_HsDm_XY_0525_broad_peaks.bed

#K27
macs2 bdgcmp -t Sample_RR_HsDm_K27Ac_0525_scalled_pileup.bdg -c Sample_RR_HsDm_K27Ac_0525_scalled_control_lambda.bdg -m qpois -o Sample_RR_HsDm_K27Ac_0525_qvalue.bed
macs2 bdgbroadcall -i Sample_RR_HsDm_K27Ac_0525_qvalue.bed -c 1.302 -C 1 -l 189 -g 25 -o  Sample_RR_HsDm_K27Ac_0525_broad_peaks.bed

#'Dux4
macs2 bdgcmp -t Sample_RR_HsDm_DUX4_0525_scalled_pileup.bdg -c Sample_RR_HsDm_DUX4_0525_scalled_control_lambda.bdg -m qpois -o Sample_RR_HsDm_DUX4_0525_qvalue.bed

macs2 bdgpeakcall -i Sample_RR_HsDm_DUX4_0525_qvalue.bdg -l 189 -g 25 -c 2 -o Sample_RR_HsDm_DUX4_0525_narrowPeaks.bed
