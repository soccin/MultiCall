#!/bin/bash

VEPROOT=/rtsess01/compute/juno/bic/ROOT/work/socci/Users/SawyersC/RomeroR1/NextFlow/MultiCall/resources/vep
perl ../../resources/vcf2maf-1.6.21/vcf2maf.pl \
    --vep-path $VEPROOT/bin \
    --vep-data $VEPROOT/vep_data \
    --vep-forks 24  \
    --ref-fasta $VEPROOT/vep_data/mus_musculus/102_GRCm38/Mus_musculus.GRCm38.dna.toplevel.fa.gz \
    --species mus_musculus \
    --ncbi-build GRCm38 \
    --input-vcf concat.vcf \
    --output-maf vep.maf \
    --tumor-id RR6_RR6 \
    --retain-info CALLER,DP \
    --retain-fmt AD
