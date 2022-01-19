#!/usr/bin/python3

import argparse
import re, sys, time

usage= '''  Description: Parse hmmsearch output '''

parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-i", "--input", required=True, help="hmmsearch domtblout format")
parser.add_argument("-v","--version", action="version", version='%(prog)s 1.4.0.21 ')
args = parser.parse_args()
parser.parse_args()

def isfloat(item):
    try:
        inS=float(item)
        return 'Yes'
    except ValueError:
        return 'No'


with open (args.input, 'r') as infile:
    for line in infile:
        if line[0]!='#':
            Line=line.rstrip()
            newLine = re.sub('\s+', '\t', Line).split('\t')
            #print(newLine[21:23])
            if newLine:
                if len(newLine)>22:
                    if isfloat(newLine[12])=='Yes':
                        description= newLine[22:]
                        desJoin= ' '.join(description)
                        print(args.input.split('/')[1][:-4], newLine[0] , newLine[3], newLine[6], newLine[12], newLine[13], newLine[17], newLine[18], newLine[19], newLine[20], desJoin, sep='\t')
                else:
                    print(args.input.split('/')[1][:-4],Line)
                    sys.exit()
            else:
                print(args.input.split('/')[1][:-4],Line)
                sys.exit()
