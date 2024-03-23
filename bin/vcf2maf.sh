#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

MROOT=$SDIR/../resources/
VEPROOT=$MROOT/vep
VCF2MAF=$MROOT/vcf2maf-1.6.21/vcf2maf.pl

IVCF=$1
NID=$2
TID=$3

perl $VCF2MAF \
    --vep-path $VEPROOT/bin \
    --vep-data $VEPROOT/vep_data \
    --vep-forks 12  \
    --ref-fasta $VEPROOT/vep_data/mus_musculus/102_GRCm38/Mus_musculus.GRCm38.dna.toplevel.fa.gz \
    --species mus_musculus \
    --ncbi-build GRCm38 \
    --retain-info CALLER,DP \
    --retain-fmt AD \
    --input-vcf $IVCF \
    --output-maf ${IVCF/.vcf/.vep.maf} \
    --normal-id $NID \
    --tumor-id $TID
