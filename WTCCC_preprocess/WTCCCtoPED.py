#!/usr/bin/python2.7

import numpy as np
import gzip
import os
import io
import sys
import cStringIO
import subprocess

def main():
    ##############################################
    #The WTCCC genotype data
    io_method = cStringIO.StringIO
    p = subprocess.Popen(["zcat",sys.argv[1]], stdout = subprocess.PIPE)
    in_file = io_method(p.communicate()[0])
    chrom = sys.argv[2]
    samp_file = open(sys.argv[3],'r')
    #Affy to RSID Conversion file
    RS_file = io.open(sys.argv[4],'r', encoding='utf-16-le')
    #List of UNIQ SNPS, pre made with BASH, should be in disease directory
    SNP_file = open('fixed'+str(chrom)+'.snps','r')


    ##############################################
    #1=control 2=case
    status = {'1958BC':1,'NBS':1,'BD':2,'CAD':2,'HT':2,'IBD':2,'RA':2,'T1D':2,'T2D':2}
    case_control = status[sys.argv[5]]
  
    #READ stuff in
    samples = np.loadtxt(samp_file,dtype='a16')
    samp_file.close()
    IDS = np.loadtxt(RS_file,skiprows=1, dtype=str)
    RS_file.close()
    SNPS = np.loadtxt(SNP_file,dtype=str)
    SNP_file.close()

    ##############################################
    #Make a list of data types for the genotype data
    #Genotype data will all be three characters -> 'a3'
    a = np.zeros((1,len(SNPS)),dtype='a3')
    a[:]='a3'

    #Types for SNP columns is a list of tuples:
    #(SNP_ID,'a3')
    typ =zip(SNPS,a[0])

    #Types of standard columns in .PED:
    stand = [('fam','i1'),('ind','a16'),('pat','i1'),('mat','i1'),('sex','i1'),('pheno','i1')]
    dash = [i for i,v in enumerate(typ) if v[0]=='---']
    counter = 0
    for i in dash:
        typ[i]=('rsbad'+str(counter),'a3')
        counter = counter +1

    #Make the array to be filled with data
    #This will, by default, fill in 0 for all numerical data type elements
    #and fill in a '' for all character type elements
    dat = np.zeros((len(samples),),dtype=stand + typ)
    
    ##############################################
    #Assign the column 'ind' to the names of samples
    dat['ind']=samples
    
    #Assign the column 'pheno' to the case_control status:
    dat['pheno']=case_control
    
    #assign all sex to female since we are not using sex at all
    dat['sex']=1
    #All other columns are currently not in need of being assigned 

    

    ##############################################
    #Make a dictionary of AFFY -> RSID
    IDict = dict(zip(IDS[:,0],IDS[:,1]))
    #Make a dictionary of WTCCID -> Position in Array
    INDict = dict(zip(samples,range(len(samples))))
     
    #For each line in the genotype data
    in_file.seek(0)
    for line in in_file:
        x = line.split()
        #If the WTCCC ID is in the Sample File
        if x[1] in INDict:
            #Then assign data to that individual:
            #If the SNP_ID needs to be translated from Affy->RSID
            if x[0] in IDict:
                if IDict[x[0]]=='---':
                    continue 
                else:
                #Then use the dictionary:
                    dat[IDict[x[0]]][INDict[x[1]]]=x[2][0]+' '+x[2][1]     
                
            else:
                #Just use the given ID
                dat[x[0]][INDict[x[1]]]=x[2][0]+' '+x[2][1]

    in_file.close()
    ##############################################
    #Fill in unassigned Genotypes with zero, which is
    #the default entry for 'missing' in PLINK! 1.07
    for i in range( len(dat) ):
        for j in range( len(dat[i]) ):
            if dat[i][j] == '':
                dat[i][j]='0 0'

    outfile = open(sys.argv[6], 'w')
    np.savetxt(outfile,dat,delimiter=' ', fmt="%s")
    outfile.close()
if __name__ == '__main__':
  main()    





 
