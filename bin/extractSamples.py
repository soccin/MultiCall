#!/usr/bin/env python3

import sys
import csv

def main(inputcsv,sampleType):
    samples=set()
    status="1" if sampleType=="TUMOR" else "0"
    with open(inputcsv) as fp:
        cin=csv.DictReader(fp,dialect="excel")
        for rec in cin:
            sid=rec["patient"],rec["sample"],rec["status"]
            if sampleType=="ALL" or status==rec["status"]:
                samples.add(sid)
            
    return samples


if __name__=="__main__":
    argv=sys.argv
    if(len(argv)<2):
        print("\n\n   usage: extractSamples.py INPUT_CSV [NORMAL|TUMOR]\n\n")
        sys.exit()
    inputcsv=argv[1]
    if(len(argv)>2):
        sampleType=argv[2].upper()
    else:
        sampleType="ALL"

    samples=main(inputcsv,sampleType)

    for si in samples:
        print(";".join([si[0],si[1]]))
