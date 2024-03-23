#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

SAMPLE=$1
VDIR=$(realpath $2)

TDIR=int/$SAMPLE
mkdir -p $TDIR

for path in $VDIR/*; do
    algo=$(basename $path)
    echo $algo
    VCF=$(find $path -name '*vcf.gz' | egrep "mutect2.filtered.vcf.gz|strelka.variants.vcf.gz|.freebayes.vcf.gz")
    $SDIR/normalizeAndTag.sh $algo $VCF | bgzip -c - > $TDIR/${algo}.vcf.gz
    tabix $TDIR/${algo}.vcf.gz
done


