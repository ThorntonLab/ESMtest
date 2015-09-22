#!/bin/sh                                                                                                                                                         
#RUNS ESM TEST; user must edit the number of CORES
#MODIFY FOR YOUR SYSTEM
#$ -q krt,bio,pub64                                                                                                                                               
#$ -N ESMTESTBIG
#$ -hold_jid PERMBIGDATA
#$ -pe openmp 4

#
#Modify as needed for your system                                                                                                                                 
module purge
module load krthornt/thorntonlab/1.0
#NEED:
#zlib
#boost
#hdf5


cd $SGE_O_WORKDIR

#Process permutations in chunks of 50 records at a time                                                                                                           
esmk -o bigfake.esmpv.txt -w 10000 -j 1000 -k 50 -n $CORES -r 0.5 --cmarkers 50 --cperms 10000 --nperms 100000 bigfake.1.perms.h5 bigfake.2.perms.h5
