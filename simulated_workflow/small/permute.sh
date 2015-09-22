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
#plink

mkfifo fake.1.mperm.dump.all
mkfifo fake.2.mperm.dump.all

#Do perms in chunks
plink --noweb --bfile fake --r2 --assoc mperm=10000 --mperm-save-all --map3 --out fake.1 --seed 1 --threads 1 &
plink --noweb --bfile fake --r2 --assoc mperm=10000 --mperm-save-all --map3 --out fake.2 --seed 2 --threads 1 &

#Process permutations in chunks of 50 records at a time
perms2h5 -i fake.1.mperm.dump.all -o fake.1.perms.h5 -b fake.bim -n 50 -l fake.1.ld
perms2h5 -i fake.2.mperm.dump.all -o fake.2.perms.h5 -b fake.bim -n 50 -l fake.1.ld
#delete named pipe
rm -f fake.*.mperm.dump.all

#Delete other needless output
rm -f fake.*.assoc.mperm

