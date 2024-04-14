#!/bin/bash

module load samtools/1.19.2

SDIR="$( cd "$( dirname "$0" )" && pwd )"

MROOT=$(realpath $SDIR/../resources)
VEPROOT=$MROOT/vep
VEPDATA=$MROOT/vep/vep_data

VCF2MAF=$MROOT/vcf2maf-1.6.21/vcf2maf.pl

IVCF=$1
NID=$2
TID=$3

perl $VCF2MAF \
    --vep-path $SDIR \
    --vep-data $VEPDATA \
    --vep-forks 12  \
    --ref-fasta $VEPDATA/mus_musculus/102_GRCm38/Mus_musculus.GRCm38.dna.toplevel.fa.gz \
    --species mus_musculus \
    --ncbi-build GRCm38 \
    --retain-info CALLER,DP,STATUS \
    --retain-fmt AD \
    --input-vcf $IVCF \
    --output-maf ${IVCF/.vcf/.vep.maf} \
    --normal-id $NID \
    --tumor-id $TID
