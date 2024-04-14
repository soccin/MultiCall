#!/bin/bash

set -eu

TAG=$1
TMP=$(mktemp -p .)
cat /dev/stdin >$TMP

cat $TMP | egrep "^##"
echo '##INFO=<ID=CALLER,Number=1,Type=String,Description="Name of mutation caller">'
cat $TMP | egrep "^#CHROM"
cat $TMP | egrep -v "^#" | awk -v tag=$TAG 'BEGIN{OFS="\t"}{$8=$8";CALLER="tag;print $0}'

rm $TMP

