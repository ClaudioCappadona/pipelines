#!/bin/bash -l

module load conda/anaconda3
source activate gatk4

wdir=/mnt/isilon/duga/joint_call_tmp/vcf/vqsr

outdir=/mnt/isilon/duga/joint_call_tmp/metrics

#vcftools --gzvcf $wdir/allChr_95_95recalibrated_demultiplex_max_mc.vcf.gz --depth --out $outdir/allChr_95_95recalibrated_demultiplex_max_mc &

#vcftools --gzvcf $wdir/allChr_95_95recalibrated_demultiplex_max_mc.vcf.gz --site-depth --out $outdir/allChr_95_95recalibrated_demultiplex_max_mc &

vcftools --gzvcf $wdir/allChr_95_95recalibrated_demultiplex_max_mc.vcf.gz --site-mean-depth --out $outdir/allChr_95_95recalibrated_demultiplex_max_mc

#vcftools --gzvcf $wdir/allChr_95_95recalibrated_demultiplex_max_mc.vcf.gz --TsTv-summary --out $outdir/allChr_95_95recalibrated_demultiplex_max_mc


exit
