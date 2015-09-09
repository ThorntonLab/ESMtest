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
    #dbSNP poitional data and Merger Archive
    pos_file = gzip.open(sys.argv[5],'r')
    merge_file = gzip.open(sys.argv[6],'r')
    disease=sys.argv[7]
    ##############################################
    #1=control 2=case
    status = {'1958BC':1,'NBS':1,'BD':2,'CAD':2,'HT':2,'IBD':2,'RA':2,'T1D':2,'T2D':2}
    case_control = status[disease]
  
    #READ stuff in
    samples = np.loadtxt(samp_file,dtype='a16')
    samp_file.close()
    IDS = np.loadtxt(RS_file,skiprows=1, dtype=str)
    RS_file.close()
    SNPS = np.loadtxt(SNP_file,dtype=str)
    SNP_file.close()

    mergeArch = {'old':'new'}

    for line in merge_file:
        x = line.split()
        mergeArch[x[0]] = x[1]

    merge_file.close()

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
    ped = np.zeros((len(samples),),dtype=stand + typ)
    
    ##############################################
    #Assign the column 'ind' to the names of samples
    ped['ind']=samples
    
    #Assign the column 'pheno' to the case_control status:
    ped['pheno']=case_control
    
    #assign all sex to female since we are not using sex at all
    ped['sex']=1
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
        #if the quality score is above 0.9
        if float(x[3]) > 0.90:
        #If the WTCCC ID is in the Sample File
            if x[1] in INDict:
            #Then assign data to that individual:
            #If the SNP_ID needs to be translated from Affy->RSID
                if x[0] in IDict:
                    if IDict[x[0]]=='---':
                        continue 
                    else:
                #Then use the dictionary:
                        ped[IDict[x[0]]][INDict[x[1]]]=x[2][0]+' '+x[2][1]     
                else:
                #Just use the given ID
                    ped[x[0]][INDict[x[1]]]=x[2][0]+' '+x[2][1]
        else:
            print >> sys.stderr, x[0:2]
    in_file.close()
    
    ##############################################
    #Fill in unassigned Genotypes with zero, which is
    #the default entry for 'missing' in PLINK! 1.07
    for i in range( len(ped) ):
        for j in range( len(ped[i]) ):
            if ped[i][j] == '':
                ped[i][j]='0 0'
    #Begin making MAP, but update 'dat' to match any SNP mergers
    merged=[]
    badmerge=[]
    nSNPs=[]
    SNP_q = {'id':'index'}
    i = 0
    #list dict's for dealing with messy mergers
    names = list(ped.dtype.names)
    flip = {'A A':'T T','A T':'T A','A G':'T C','A C':'T G','T T':'A A','T A':'A T','T G':'A C','T C':'A G','G G':'C C','G A':'C T','G T':'C A','G C':'C G','C C':'G G','C A':'G T','C T':'G A','C G':'G C','0 0':'0 0'}
    #make dictionary of rsID -> position in MAP                                 
    for a in SNPS:
        #If SNP has been merged
        if a.strip('rs')in mergeArch:
            b= mergeArch[a.strip('rs')]
            rsb = 'rs' + b
            #is the SNP that SNP a was merged to( SNP b ) in the original list 
            if rsb in SNPS:
                new = ped[rsb]
                old = ped[a]
                #If they are completely equal, this is good
                if set(old)==set(new):
                    names.remove(a)
                    merged.append(a + ' ' + rsb + '\n')
                else:
                    #Check for strand flip
                    check= [flip[old[j]]==new[j] for j in range(len(old))]
                    #If it is an exact match strand flip...this is very improbable, given that the SNPs were merged
                    if sum(check)==len(old):
                        names.remove(a)
                        merged.append(a + ' ' + rsb + '\n')
                    else:
                        names.remove(a)
                        names.remove(rsb)
                        badmerge.append(rsb)
                        #if we already added b to the list
                        if b in SNP_q:
                            #remove b from the list
                            nSNPs.remove(rsb)
                            #and reupdate the indexing
                            pos_b = SNP_q[b]
                            del SNP_q[b]
                            for pos_snp in range(pos_b,len(SNP_q)):
                                SNP_q[SNP_q.keys()[pos_snp]]-=1
            else:
                nSNPs.append(rsb)
                SNP_q[b]=i
                i = i + 1
                merged.append(a + ' ' + rsb + '\n')
        else:
            if a not in badmerge:
                nSNPs.append(a)
                SNP_q[a.strip('rs')]=i
                i = i + 1
    
    ped = ped[names]
    MAP = np.zeros((len(nSNPs),),dtype=[('chrom','a2'),('rsid','a16'),('pos','a16')])
    MAP['chrom']=str(chrom)
    MAP['rsid']=nSNPs
    exclude = []
    for line in pos_file:
        x = line.split()
        if x[0] in SNP_q:
            if len(x)<3:
                MAP[SNP_q[x[0]]][2]= '1'
                exclude.append('rs'+x[0]+'\n')
            else:
                MAP[SNP_q[x[0]]][2]=x[2]
                
    i = 0
    for line in MAP:
        if len(line[2])==0:
            MAP[i][2]='1'
        i = i +1
   

    mapfile = open(disease+chrom+".map", 'w')
    np.savetxt(mapfile,MAP,delimiter=' ', fmt="%s")
    mapfile.close()

    excludefile = open(chrom+'exclusion.txt','w')
    excludefile.writelines(exclude)
    excludefile.close()
    
    mergefile = open(disease+chrom+'merged.txt','w')
    mergefile.writelines(merged)
    mergefile.close()

    badmergefile = open(disease+chrom+'badmerge.txt','w')
    badmergefile.writelines(badmerge)
    badmergefile.close()

    pedfile = open(disease+chrom+".ped", 'w')
    np.savetxt(pedfile,ped,delimiter=' ', fmt="%s")
    pedfile.close()
if __name__ == '__main__':
  main()    





 
