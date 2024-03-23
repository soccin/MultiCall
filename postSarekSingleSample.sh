#!/bin/bash

set -eu

SDIR="$( cd "$( dirname "$0" )" && pwd )"
export PATH=$SDIR/bin:$PATH

TARGET_BED=$SDIR/targets/M-IMPACT_v2_mm10_targets_plus15bp.bed

if [ "$#" != "2" ]; then
    echo -e "\n\n   usage: postSarekSingleSample.sh OUTPUT_DIRECTORY INPUT_CSV\n\n"
    exit
fi

IDIR=$1
INPUTCSV=$2

TUMOR=$(extractSamples.py $INPUTCSV tumor)
NORMAL=$(extractSamples.py $INPUTCSV normal)

PID=${NORMAL%%;*}
TID=${TUMOR##$PID;}
NID=${NORMAL##$PID;}

exec >> ${PID}_${NID}.log
exec 2>&1

echo PID=$PID
echo TID=$TID 
echo NID=$NID

NTAG=${NORMAL/;/_}
TTAG=${TUMOR/;/_}

VCDIR=$(find $IDIR -type d | egrep "/${PID}/variant_calling$")

MUTECT_VCF=$VCDIR/mutect2/${TID}_vs_${NID}/*.mutect2.filtered.vcf.gz
FREEBAYES_VCF=$VCDIR/freebayes/${TID}_vs_${NID}/*.freebayes.vcf.gz
STRELKA_SNVS_VCF=$VCDIR/strelka/${TID}_vs_${NID}/*.somatic_snvs.vcf.gz
STRELKA_INDELS_VCF=$VCDIR/strelka/${TID}_vs_${NID}/*.somatic_indels.vcf.gz

#
# Caller specific post process to
# - fix sample names (NORMAL/TUMOR in strelka)
# - merge snps/indels (fu strelka)
# - fix another caller specific nonsense
#

ODIR=int/${TID}___${NID}
mkdir -p $ODIR

post_strelka.sh \
    $STRELKA_SNVS_VCF $STRELKA_INDELS_VCF \
    $NTAG $TTAG $TARGET_BED \
    > $ODIR/strelka.vcf
vcf2maf.sh $ODIR/strelka.vcf $NTAG $TTAG &

normalizeAndTag.sh mutect2 $MUTECT_VCF >$ODIR/mutect2.vcf
vcf2maf.sh $ODIR/mutect2.vcf $NTAG $TTAG &

#
# For Freebayes remove samples with null calls ".:.:.:."
#
normalizeAndTag.sh freebayes $FREEBAYES_VCF \
    | fgrep -v ".:.:." \
    >$ODIR/freebayes.vcf
vcf2maf.sh $ODIR/freebayes.vcf $NTAG $TTAG

wait

Rscript mergeMAFs.R $ODIR/*.vep.maf

echo DONE $INPUTCSV
