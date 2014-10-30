#!/bin/sh

#$ -N MAKEDATA
#$ -q krt,bio,pub64
#uses ms to simulate a data set of 500 SNPs in 6,000 individuals


#Modify as needed for your GE system
module load krthornt/libsequence/1.8.0
module load zlib/1.2.7
module load plink/1.90a

SEED1=$RANDOM
SEED2=$RANDOM
SEED3=$RANDOM

#Simulate the data with a recent severe bottleneck, for no reason other than that we can
~/msdir/ms 6000 1 -s 500 -eN 0.001 0.1 -eN 0.01 1 -seed $SEED1 $SEED2 $SEED3 | ~/ESMtest/fake_data/ms2plink fake.ped fake.map

#Make the binary input files for plink.  Apply HWE filters, exclusions of SNPS, etc., at this stage
#We apply liberal HWE filter here just so that something does get excluded.
plink --noweb --file fake --make-bed --map3 --out fake --silent --hwe 0.1

#The .ped file is no longer needed
rm -f fake.ped

#delete the log file
rm -f fake.log