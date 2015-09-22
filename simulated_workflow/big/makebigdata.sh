#!/bin/sh
#MODIFY FOR YOUR CLUSTER!
#IF YOU ARE RUNNING THIS ON A SINGLE MACHINE THEN
#REMOVE THIS
#$ -N MAKEBIGDATA
#$ -q krt,bio,pub64
#uses ms to simulate a data set of 40000 SNPs in 6,000 individuals
########

#Modify as needed for your GE system
module purge
module load krthornt/thorntonlab/1.0
module load plink/1.90a
#NEED:
#zlib
#boost
#hdf5
#plink

SEED1=$RANDOM
SEED2=$RANDOM
SEED3=$RANDOM

#Simulate the data with a recent severe bottleneck, for no reason other than that we can
#Simulate 40K SNPS
~/apps/msdir/ms 6000 1 -s 40000 -eN 0.001 0.1 -eN 0.01 1 -seed $SEED1 $SEED2 $SEED3 | ms2plink bigfake.ped bigfake.map

#Make the binary input files for plink.  Apply HWE filters, exclusions of SNPS, etc., at this stage
#We apply liberal HWE filter here just so that something does get excluded.
plink --noweb --file bigfake --make-bed --map3 --out bigfake --silent --hwe 0.1
