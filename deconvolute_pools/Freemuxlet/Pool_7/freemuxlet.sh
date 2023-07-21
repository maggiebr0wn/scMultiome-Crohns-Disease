#!/usr/bin/bash

BAMFILE="/storage/home/mfisher42/scProjects/CD_Subra/rawdata/Maggie_Organized_Data/cellranger_output/Pool_7/outs/gex_possorted_bam.bam"
BARCODES="/storage/home/mfisher42/scProjects/CD_Subra/rawdata/Maggie_Organized_Data/cellranger_output/Pool_7/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"

# 1.) need to run dsc-pileup, a software tool to pileup reads and base quality for each overlapping SNPs and each barcode. Allows us to run demuxlet/freemuxlet without going over the bamfile again. Requires: SAM/BAM/CRAM file produced from 10x, a group list file to specify valid droplet set, a VCF file containing AC and AN from reference panel (ie. 1000g)

# This step took 14 hours; maybe try to do each chromosome in parallel next time

mkdir tmp
../sort_vcf_same_as_bam.sh $BAMFLE ../ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz > ./tmp/samples.sorted_as_in_bam.vcf

popscle dsc-pileup --sam $BAMFILE --group-list $BARCODES --vcf tmp/samples.sorted_as_in_bam.vcf --out samples_to_demultiplex.pileup

# 2.) Freemuxlet is then run in cases where genotype data is unavailable; used to deconvolute sample identity and identify multiplets when multiple samples are pooled by barcoded single cell sequencing.

# This step took ~30% Memory
popscle freemuxlet --plp samples_to_demultiplex.pileup --out freemuxlet_output --nsample 3
