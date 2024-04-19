#!/bin/bash

set -eu

module load bcftools/1.19

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

TMP=$(mktemp -p .)
trap "rm ${TMP}*" EXIT

#
# bcftools view -R BED needs an indexed vcf
# bedtools intersect does not work with all vcf's
#

bcftools concat $SNVS $INDELS -a | bgzip -c - >${TMP}.vcf.gz
tabix -p vcf ${TMP}.vcf.gz

bcftools view -R $TARGET_BED ${TMP}.vcf.gz \
    | bcftools sort - \
    | bcftools norm -m- \
    | perl -pe 's/NORMAL/'${NORMAL}'/ if /^#C/; s/TUMOR/'${TUMOR}'/ if /^#C/' \
    | tagVCF.sh strelka
