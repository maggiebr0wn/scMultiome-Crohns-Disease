#!/usr/bin/bash

plink_vcf="/storage/home/mfisher42/scProjects/CD_Subra/rawdata/Genotype_Data/Genos_by_Pool/Pool7/LiftOver_to_Hg38/hg38_output.vcf"

# in dir with file, fix order of vcf chromosome to be lexicographically
cat $plink_vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > sorted_vcf.vcf

# run Demuxlet
demuxlet --sam sorted_bam.bam --tag-group CB --tag-UMI UB --vcf sorted_vcf.vcf --field GT --out pool7_demux
