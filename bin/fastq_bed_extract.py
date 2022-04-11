#!/usr/bin/env python3
__author__ = "Victor de Jager"
__copyright__ = "Copyright 2022, WUR-MIB"
__credits__ = ["Victor de Jager"]
__license__ = "CC0"
__version__ = "1.0.0"
__status__ = "Development"

import os,sys,re
import string
import argparse

def bed_to_dict(bed):
    beddict={}
    with open(bed,'r') as fh:
        for line in fh:
            line = line.strip().split()
            beddict["@%s"%line[0]]={'start':line[1],
                          'end':line[2],
                          'name':line[3],
                          'qual':line[4],
                          'strand':line[5]}

    return beddict

def reverse_complement_table(seq):
    matrix = str.maketrans("ACTG", "TGAC")
    return seq.translate(matrix)[::-1]


def process(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    d = {k: v for k, v in zip(ks, lines)}
    d['name'] = d['name'].split()[0]
    return d


def process_fastq_bed(fastq, bed, out):

    #make dictionary from bedfile
    beddict = bed_to_dict(bed)

    # empty list
    filtered_records= []
    n = 4
    l=0

    with open(fastq, 'r') as fh:
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            if len(lines) == n:
                record = process(lines)
                l=l+1

                if record['name'] in beddict:
                    b = beddict[ record['name'] ]
                    start = int(b['start'])
                    end = int(b['end'])
                    name = b['name']
                    strand = b['strand']
                    recordname = record['name']
                    
                    if (strand=="-"):
                        record['sequence'] = reverse_complement_table(record['sequence'][start:end])
                    else:
                        record['sequence'] = record['sequence'][start:end]
                    
                    record['quality'] = record['quality'][start:end]
                    record['name'] = f"@{name}::{recordname}:{start}:{end}({strand})"
                    
                    filtered_records.append(record)
                lines = []
    with open(out,'w') as ofh:
        for fr in filtered_records:
            ofh.write(fr['name'])
            ofh.write("\n")
            ofh.write(fr['sequence'])
            ofh.write("\n")
            ofh.write(fr['optional'])
            ofh.write("\n")
            ofh.write(fr['quality'])
            ofh.write("\n")
        ofh.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq", help="Fastq file to be filtered", required=True)
    parser.add_argument("--bed", help="bedfile to use for filtering", required=True)
    parser.add_argument("--out", help="Fastq output file", required=True)
   

    args = parser.parse_args()
    if not args.fastq and not args.bed and not args.out:
        print("All arguments are required")
        sys.exit(1)

    # Process the fastq and bed file
    fastq = args.fastq
    bed = args.bed
    out = args.out
    process_fastq_bed(fastq, bed, out)


