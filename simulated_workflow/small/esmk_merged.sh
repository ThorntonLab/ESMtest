#!/bin/sh                                                                                                                                                         

#Modify as needed for your system                                                                                                                                 
module purge
module load krthornt/thorntonlab/1.0
module load plink/1.90a
#NEED:
#zlib
#boost
#hdf5
#plink


#Process permutations in chunks of 50 records at a time                                                                                                           
esmk -o fake_merged.esmpv.txt -w 10000 -j 1000 -k 50 -n 1 -r 0.5 --cmarkers 50 --cperms 10000 --nperms 10000 fake.all.perms.h5
