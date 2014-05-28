#!/bin/sh

#Does 500 perms of the data 5 times.
#This mimics lots of perms done in chunks

#$ -N PERMDATA
#$ -hold_jid MAKEDATA
#$ -t 1-5

#Modify as needed for your system
module load plink/1.90a
module load zlib/1.2.7
module load boost/1.53.0

cd $SGE_O_WORKDIR

#We use a named pipe to skip writing the permutation output in plain text
mkfifo fake.$SGE_TASK_ID.mperm.dump.all

SEED=$SGE_TASK_ID
#Do this chunk of perms
plink --noweb --bfile fake --assoc mperm=500 --map3 --mperm-save-all --hwe 1e-6 --out fake.$SGE_TASK_ID  --seed $SEED &

$1/processperms/perms2gz -i fake.$SGE_TASK_ID.mperm.dump.all -o fake.$SGE_TASK_ID.mperm.dump.all.gz 

#delete named pipe
rm -f fake.$SGE_TASK_ID.mperm.dump.all

