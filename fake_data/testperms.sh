#!/bin/bash

#$ -q krt
#$ -t 1-64

module load krthornt/plink/1.07
module load krthornt/libsequence/1.8.0

cd $SGE_O_WORKDIR
SEED1=`echo "$SGE_TASK_ID*$RANDOM"|bc -l`
SEED3=`echo "$SGE_TASK_ID*$RANDOM"|bc -l`
SEED2=`echo "$SGE_TASK_ID*$RANDOM"|bc -l`

ms 12000 1 -s 250 -seed $SEED1 $SEED2 $SEED3 | ~/src/ESMtest/fake_data/ms2plink testdata.$SGE_TASK_ID.ped testdata.$SGE_TASK_ID.map
#Make binary input files
plink --noweb --file testdata.$SGE_TASK_ID --make-bed --map3 --out testdata.$SGE_TASK_ID --silent
#Do perms
/usr/bin/time -f "%e %M" -o ptime1million.$SGE_TASK_ID.txt plink --noweb --bfile testdata.$SGE_TASK_ID --assoc --map3 --mperm 1000000 --mperm-save-all --hwe 1e-6 --out testdata.$SGE_TASK_ID --silent

#We cannot suppress the log file, but we can delete it:
rm -f testdata.$SGE_TASK_ID.log

#After running a job like this, the original .bed and .map files should probably be gzipped and tarred up

