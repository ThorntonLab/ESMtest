#!/usr/bin/python2.7

import numpy as np
import gzip
import os
import io
import sys

def main():
    #The WTCCC genotype data
    pos_file = gzip.open(sys.argv[1],'r')
    merge_file = gzip.open(sys.argv[2],'r')
    chrom = sys.argv[3]

    #List of UNIQ SNPS, pre made with BASH, should be in disease directory
    SNP_file = open('fixed'+str(chrom)+'.snps','r')

    #READ stuff in
    mergeArch = {'old':'new'}

    for line in merge_file:
        x = line.split()
        mergeArch[x[0]] = x[1]

    merge_file.close()

    SNPS= SNP_file.read().splitlines()    
    SNP_file.close()

    SNP_q = {'id':'index'}
    i = 0
    for a in SNPS:
        if a.strip('rs')in mergeArch:
            b= mergeArch[a.strip('rs')]
            SNP_q[b]=i
        else:
            SNP_q[a.strip('rs')]=i
        i = i + 1

    dat = np.zeros((len(SNPS),),dtype=[('chrom','a2'),('rsid','a16'),('pos','a16')])
    dat['chrom']=str(chrom)
    dat['rsid']=SNPS
    exclude = []
    for line in pos_file:
        x = line.split()
        if x[0] in SNP_q:
            if len(x)<3:
                dat[SNP_q[x[0]]][2]= '1'
                exclude.append('rs'+x[0]+'\n')
            else:
                dat[SNP_q[x[0]]][2]=x[2]



    outfile = open(sys.argv[4], 'w')

    np.savetxt(outfile,dat,delimiter=' ', fmt="%s")

    outfile.close()
    excludefile = open(chrom+'exclusion.txt','w')
    excludefile.writelines(exclude)
    excludefile.close()

if __name__ == '__main__':
  main()






 
