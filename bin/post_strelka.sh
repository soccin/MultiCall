#!/bin/bash

SNVS=$1
INDELS=$2
NORMAL=$3
TUMOR=$4
TARGET_BED=$5

#
# Strelka is intentionally difficult
# - separates snvs from indels
# - mislabels samples
#

bcftools concat $SNVS $INDELS -a \
    | bedtools intersect -a - -b $TARGET_BED -wa -header \
    | bcftools sort - \
    | bcftools norm -m- \
    | perl -pe 's/NORMAL/'${NORMAL}'/ if /^#C/; s/TUMOR/'${TUMOR}'/ if /^#C/' \
    | tagVCF.sh strelka


