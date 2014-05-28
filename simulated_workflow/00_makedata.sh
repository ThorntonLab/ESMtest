#!/bin/sh

#$ -N MAKEDATA

#uses ms to simulate a data set of 500 SNPs in 6,000 individuals

cd $SGE_O_WORKDIR

#Modify as needed for your GE system
module load krthornt/libsequence/1.8.0
module load zlib/1.2.7
module load plink/1.90a

SEED1=$RANDOM
SEED2=$RANDOM
SEED3=$RANDOM

#Simulate the data
ms 6000 1 -s 500 -seed $SEED1 $SEED2 $SEED3 | $1/fake_data/ms2plink fake.ped fake.map

#Make the binary input files for plink
plink --noweb --file fake --make-bed --map3 --out fake --silent

#The .ped file is no longer needed
rm -f fake.ped