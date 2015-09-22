#!/bin/sh

module purge
module load krthornt/thorntonlab/1.0
module load plink/1.90a
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

