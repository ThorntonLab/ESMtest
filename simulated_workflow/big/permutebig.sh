#!/bin/sh

#Does 5 perms of the data 5 times.
#This mimics lots of perms done in chunks
##MODIFY FOR YOUR CLUSER
#IF YOU ARE RUNNING ON A REGULAR MACHINE
#THEN REMOVE THIS
#$ -q krt,bio,pub64,free64
#$ -N PERMBIGDATA
#$ -hold_jid MAKEBIGDATA
#$ -t 1-2
#$ -o /dev/null
#$ -pe openmp 2
#Modify as needed for your system
module purge
module load krthornt/thorntonlab/1.0
module load plink/1.90a
#NEED:
#zlib
#boost
#hdf5
#plink

cd $SGE_O_WORKDIR

mkfifo bigfake.$SGE_TASK_ID.mperm.dump.all

SEED=$SGE_TASK_ID
#Do this chunk of perms
plink --noweb --bfile bigfake --r2 --assoc mperm=50000 --mperm-save-all --map3 --out bigfake.$SGE_TASK_ID  --seed $SEED --threads $CORES &
#Process permutations in chunks of 50 records at a time
perms2h5 -i bigfake.$SGE_TASK_ID.mperm.dump.all -o bigfake.$SGE_TASK_ID.perms.h5 -b bigfake.bim -n 50 -l bigfake.${SGE_TASK_ID}.ld
#delete named pipe
rm -f bigfake.$SGE_TASK_ID.mperm.dump.all

#Delete other needless output
rm -f bigfake.$SGE_TASK_ID.assoc.mperm
#rm -f bigfake.$SGE_TASK_ID.assoc
