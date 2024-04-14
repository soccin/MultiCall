#!/bin/bash

set -eu

SDIR="$( cd "$( dirname "$0" )" && pwd )"

module load bcftools/1.19

TARGETBED=$SDIR/../targets/M-IMPACT_v2_mm10_targets_plus15bp.bed

TMP=$(mktemp -p .)

TAG=$1
IVCF=$2

#
# - bedtools does not work with VARDICT VCF files (<DEL> events)
#
# - tabix returns if the same event hits multiple targets
#
# - bcftools seems to do the best job although it does mess with
# floating point numbers (11.000 => 11)
#

bcftools view -R $TARGETBED $IVCF \
    | bcftools norm -m- \
    > $TMP

cat $TMP | egrep "^##"
echo '##INFO=<ID=CALLER,Number=1,Type=String,Description="Name of mutation caller">'
cat $TMP | egrep "^#CHROM"
cat $TMP | egrep -v "^#" | awk -v tag=$TAG 'BEGIN{OFS="\t"}{$8=$8";CALLER="tag;print $0}'

rm $TMP

