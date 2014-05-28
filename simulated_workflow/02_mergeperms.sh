#!/bin/sh

#Merge all the permuted data into 1 single file

#$ -N MERGEPERMS
#$ -hold_jid PERMDATA

module load boost/1.53.0
module load zlib/1.2.7

cd $SGE_O_WORKDIR

$1/processperms/mergeperms -o fake_perms.merged.gz -m fake.map fake.*.mperm.dump.all.gz 

rm -f fake.*.mperm.dump.all.gz