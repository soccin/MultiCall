#!/bin/bash

if [ "$#" != "2" ]; then
    echo -e "\n\n   usage: postSarekSingleSample.sh OUTPUT_DIRECTORY INPUT_CSV\n\n"
    exit
fi

ODIR=$1
INPUTCSV=$2

