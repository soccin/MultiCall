#!/bin/bash

set -eu


SDIR="$( cd "$( dirname "$0" )" && pwd )"
export PATH=$SDIR/bin:$PATH

VEPVERSION=$(vep --help | fgrep ensembl-vep | awk '{print $3}')

if [ "$VEPVERSION" != "102.0" ]; then
    echo -e "\n\n   vep not properly installed; see instructions\n"
    echo -e "VEPVERSION=[${VEPVERSION}]\n\n"
    exit 1
fi


TARGET_BED=$SDIR/targets/M-IMPACT_v2_mm10_targets_plus15bp.bed

if [ "$#" != "4" ]; then
    echo -e "\n   usage: postSarekPair.sh VAR_CALL_DIR PID NID TID\n"
    echo -e "        VEPVERSION=[${VEPVERSION}]\n\n"
    exit
fi


VCDIR=$(realpath $1)
PID=$2
NID=$3
TID=$4

# exec >> ${PID}_${NID}.log
# exec 2>&1

NTAG=${PID}_${NID}
TTAG=${PID}_${TID}

MUTECT_VCF=$VCDIR/mutect2/${TID}_vs_${NID}/*.mutect2.filtered.vcf.gz
FREEBAYES_VCF=$VCDIR/freebayes/${TID}_vs_${NID}/*.freebayes.vcf.gz
STRELKA_SNVS_VCF=$VCDIR/strelka/${TID}_vs_${NID}/*.somatic_snvs.vcf.gz
STRELKA_INDELS_VCF=$VCDIR/strelka/${TID}_vs_${NID}/*.somatic_indels.vcf.gz

#
# Vardict output not in out folder as not part of sarek
# look for it in post folder
#

VARDICT_VCF=post/variant_calling/vardict/${TID}_vs_${NID}/*.vardict.vcf.gz

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
vcf2maf.sh $ODIR/freebayes.vcf $NTAG $TTAG &

#
# For Vardict get rid of <DEL> events
#
if [ -e $VARDICT_VCF ]; then
    echo Processing VARDICT VCF $(basename $VARDICT_VCF)
    normalizeAndTag.sh vardict $VARDICT_VCF \
        | fgrep -v "<DEL>" \
        >$ODIR/vardict.vcf
    vcf2maf.sh $ODIR/vardict.vcf $NTAG $TTAG
fi

wait

Rscript $SDIR/mergeMAFs.R $ODIR/*.vep.maf

echo DONE $*
