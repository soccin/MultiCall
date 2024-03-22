#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

TARGETBED=$SDIR/targets/M-IMPACT_v2_mm10_targets_plus15bp.bed

TMP=$(mktemp -p .)

TAG=$1
IVCF=$2

bedtools intersect -a $IVCF -b $TARGETBED -wa -header \
    | bcftools norm -m- \
    > $TMP

cat $TMP | egrep "^##"
echo '##INFO=<ID=CALLER,Number=1,Type=String,Description="Name of mutation caller">'
cat $TMP | egrep "^#CHROM"
cat $TMP | egrep -v "^#" | awk -v tag=$TAG 'BEGIN{OFS="\t"}{$8=$8";CALLER="tag;print $0}'

rm $TMP

