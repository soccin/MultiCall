#!/bin/bash

set -eu

export SDIR="$( cd "$( dirname "$0" )" && pwd )"
export PATH=$SDIR/bin:$PATH

VEPVERSION=$(vep --help | fgrep ensembl-vep | awk '{print $3}')

if [ "$VEPVERSION" != "102.0" ]; then
    echo -e "\n\n   vep not properly installed; see instructions\n"
    echo -e "VEPVERSION=[${VEPVERSION}]\n\n"
    exit 1
fi

if [ "$#" != "5" ]; then
    echo -e "\n   usage: postSarekPair.sh TARGETS VAR_CALL_DIR PID NID TID\n"
    echo -e "        VEPVERSION=[${VEPVERSION}]\n\n"
    exit
fi


TARGETS=$1
if [ ! -e $SDIR/targets/$TARGETS/target.resources.sh ]; then 
    echo -e "\n\nMissing target resources\n"
    echo "TARGETS=[${TARGETS}]"
    exit 1
fi

#
# This will (should) set the following variables:
# - TARGET_BED
#
source $SDIR/targets/$TARGETS/target.resources.sh

VCDIR=$(realpath $2)
PID=$3
NID=$4
TID=$5

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

VARDICT_VCF=$(find post/variant_calling/vardict -name "${TID}_vs_${NID}*.vardict.vcf.gz")

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

normalizeAndTag.sh $TARGET_BED mutect2 $MUTECT_VCF >$ODIR/mutect2.vcf
vcf2maf.sh $ODIR/mutect2.vcf $NTAG $TTAG &

#
# For Freebayes remove samples with null calls ".:.:.:."
#
normalizeAndTag.sh $TARGET_BED freebayes $FREEBAYES_VCF \
    | fgrep -v ".:.:." \
    >$ODIR/freebayes.vcf
vcf2maf.sh $ODIR/freebayes.vcf $NTAG $TTAG &

#
# For Vardict get rid of <DEL> events
#
if [ -e $VARDICT_VCF ]; then
    echo Processing VARDICT VCF $(basename $VARDICT_VCF)
    normalizeAndTag.sh $TARGET_BED vardict $VARDICT_VCF \
        | fgrep -v "<DEL>" \
        >$ODIR/vardict.vcf
    vcf2maf.sh $ODIR/vardict.vcf $NTAG $TTAG
fi

wait

Rscript $SDIR/mergeMAFs.R $ODIR/*.vep.maf

echo DONE $*
