#!/bin/sh                                                                                                                                                         
#RUNS ESM TEST; user must edit the number of CORES
#$ -q krt,bio,pub64                                                                                                                                               
#$ -N ESMTEST
#$ -hold_jid PERMDATA
#$ -pe openmp 1

#
#Modify as needed for your system                                                                                                                                 
module load zlib/1.2.7
module load boost/1.54.0
module load hdf5/1.8.11

cd $SGE_O_WORKDIR

#Process permutations in chunks of 50 records at a time                                                                                                           
/data/users/jsanjak/ESMtest/src/esmk -o fake.esmpv.txt -w 10000 -j 1000 -k 50 -n $CORES -r 0.5 --cmarkers 50 --cperms 10000 --nperms 100000 fake.1.perms.h5 fake.2.perms.h5
