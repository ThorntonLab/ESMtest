#!/usr/bin/python2.7.2
import io
import sys
import numpy as np

def main():
    
    RS_file = io.open('/bio/jsanjak/thornton/GWAS/AFFY_RSID.txt','r', encoding='utf-16-le')
    
    #List of UNIQ SNPS, pre made with BASH
    SNP_file = open(sys.argv[1],'r')
    outfile = open('fixed'+sys.argv[1],'w')
    
    IDS = np.loadtxt(RS_file,skiprows=1, dtype=str)
    RS_file.close()
    SNPS = np.loadtxt(SNP_file,dtype=str)
    SNP_file.close()

    for i in range(len(SNPS)):
        idx = np.where(IDS[:,0]==SNPS[i])[0]
        if len(idx)>0:
            outfile.write(IDS[idx[0],1]+'\n')
        else:
            outfile.write(SNPS[i]+'\n')


if __name__ == '__main__':
  main() 
