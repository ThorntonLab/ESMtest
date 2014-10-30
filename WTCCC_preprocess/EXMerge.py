#!/usr/bin/python2.7

import numpy as np
import gzip
import os
import io
import sys

def main():
    #MergeArch
    merge_file = gzip.open(sys.argv[1],'r')
    #List of UNIQ SNPS, pre made with BASH, should be in disease directory
    SNP_file = open(sys.argv[2],'r')

    #READ stuff in
    mergeArch = {'old':'new'}

    for line in merge_file:
        x = line.split()
        mergeArch[x[0]] = x[1]

    merge_file.close()

    outfile = open(sys.argv[3], 'w')
    for SNP in SNP_file:
        merged = mergeArch[SNP]
        outfile.write(merged+'\n')
    

    outfile.close()
    

if __name__ == '__main__':
  main()






 
