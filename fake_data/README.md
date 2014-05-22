#A tool to take [ms](http://home.uchicago.edu/rhudson1/source/mksamples.html) output and create a PED/MAP file pair for use with PLINK.

#Example use:  make 3,000 cases + 3,000 controls and 500 SNPs.  That is 6,000 diploids, so 12,000 chromosomes.  Using ms:

```{sh}
ms 12000 1 -s 500 | ./ms2plink fake.ped fake.map
```

##Then, to get $10^6$ permutations:

```{sh}
#We use --map3 b/c we do not have the 3rd column in the map file, which is position in cM.
plink --file fake --assoc --map3 --mperm 1000000 --mperm-save-all
```


Full example of a workflow:

```{sh}
#Simulate some data.
ms 12000 1 -s 500 | ./ms2plink fake.ped fake.map
#Make binary input files
../../plink_test/plink-1.07-srcKRT/plink --file fake --make-bed --map3 --out fake
#Do perms
../../plink_test/plink-1.07-srcKRT/plink --bfile fake --assoc --map3 --mperm 200000 --mperm-save-all --hwe 1e-6 --out fake
```

The above are executed in the script called fakeit.sh.


#Requirements on HPC:

```
module load krthornt/libsequence/1.8.0
```

#Results

I ran the script testperms.sh, which permuts 250 snps $10^6$ times.

```
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
```

The distribution of run times looks like:

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1.png) 


The mean run time is 0.9671 days, with a max of 1.0023 days.

Given that we will need to do $3 \times 10^6$ perms, and we have a 3 day run-time limit on the pub64 queue, it is probably safest to limit the permutations to 175 or 200 SNP chunkks (as run time should be linear in the number of markers).
