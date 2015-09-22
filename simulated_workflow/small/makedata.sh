#!/bin/sh

module purge
module load krthornt/thorntonlab/1.0
module load plink/1.90a
#NEED:
#zlib
#boost
#hdf5
#plink

#Simulate the data with a recent severe bottleneck, for no reason other than that we can
~/apps/msdir/ms 6000 1 -s 500 -eN 0.001 0.1 -eN 0.01 1 -seed 1 2 3 | ms2plink fake.ped fake.map

#Make the binary input files for plink.  Apply HWE filters, exclusions of SNPS, etc., at this stage
#We apply liberal HWE filter here just so that something does get excluded.
plink --noweb --file fake --make-bed --map3 --out fake --silent --hwe 0.1

