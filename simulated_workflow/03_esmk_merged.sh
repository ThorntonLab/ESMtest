#!/bin/sh                                                                                                                                                         
#RUNS ESM TEST; user must edit the number of CORES
#$ -q krt,bio,pub64                                                                                                                                               
#$ -N ESMTEST
#$ -hold_jid MERGE
#$ -pe openmp 1

#
#Modify as needed for your system                                                                                                                                 
module purge
module load krthornt/thorntonlab/1.0
module load plink/1.90a
#NEED:
#zlib
#boost
#hdf5
#plink

cd $SGE_O_WORKDIR

#Process permutations in chunks of 50 records at a time                                                                                                           
/data/users/jsanjak/ESMtest/src/esmk -o fake_merged.esmpv.txt -w 10000 -j 1000 -k 50 -n $CORES -r 0.5 --cmarkers 50 --cperms 10000 --nperms 100000 fake.all.perms.h5
