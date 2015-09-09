#!/bin/sh                                                                                                   
#$ -q krt,bio,pub64,free64
#$ -N esmtest
#$ -pe openmp 16

module load zlib/1.2.7
module load boost/1.54.0
module load hdf5/1.8.11

cd $SGE_O_WORKDIR

#Process permutations in chunks of 50 records at a time                                                                                                           
/data/users/jsanjak/ESMtest/src/esmk -o multi.16.esmpv.txt -w 10000 -j 1000 -k 50 -n 16 -r 0.5 fake.1.perms.h5 fake.2.perms.h5

