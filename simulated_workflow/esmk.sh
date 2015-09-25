#!/bin/sh                                                                                                                                         
#On UCI HPC we use                
#module purge
#module load krthornt/thorntonlab/1.0
#WHICH LOADS:
##gcc/4.8.4
##gsl/1.16
##htslib/1.2.1
##zlib/1.2.8
##boost/1.59.0
##hdf5/1.8.15-patch1
##python/2.7.10
##krthornt/libsequence/1.8.5

#module load plink/1.90a

#Modify as needed for your system                                                                                                                 
#NEED:
#zlib
#boost
#hdf5

#Process permutations in chunks of 50 records at a time                                                                                           
esmk -o fake.esmpv.txt -w 10000 -j 1000 -k 50 -n 1 -r 0.5 --cmarkers 50 --cperms 10000 --nperms 20000 fake.1.perms.h5 fake.2.perms.h5
