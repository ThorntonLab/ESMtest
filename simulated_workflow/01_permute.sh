#!/bin/sh

#Does 5 perms of the data 5 times.
#This mimics lots of perms done in chunks
#$ -q krt,bio,pub64,free64
#$ -N PERMDATA
#$ -hold_jid MAKEDATA
#$ -t 1-2
#$ -o /dev/null

#$ -pe openmp 2
#Modify as needed for your system
module load plink/1.90a
module load zlib/1.2.7
module load boost/1.54.0
module load hdf5/1.8.11

cd $SGE_O_WORKDIR

mkfifo fake.$SGE_TASK_ID.mperm.dump.all

SEED=$SGE_TASK_ID
#Do this chunk of perms
plink --noweb --bfile fake --r2 --assoc mperm=50000 --mperm-save-all --map3 --out fake.$SGE_TASK_ID  --seed $SEED --threads $CORES &
#Process permutations in chunks of 50 records at a time
~/ESMtest/src/perms2h5 -i fake.$SGE_TASK_ID.mperm.dump.all -o fake.$SGE_TASK_ID.perms.h5 -b fake.bim -n 50 -l fake.${SGE_TASK_ID}.ld
#delete named pipe
rm -f fake.$SGE_TASK_ID.mperm.dump.all

#Delete other needless output
rm -f fake.$SGE_TASK_ID.log
rm -f fake.$SGE_TASK_ID.assoc.mperm
#rm -f fake.$SGE_TASK_ID.assoc