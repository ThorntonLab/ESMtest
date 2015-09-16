#!/bin/sh
#$ -q krt,bio,pub64
#$ -N MERGE
#$ -hold_jid PERMDATA

if [ -a fake_merged.all.perms.h5 ]
then
    rm fake_merged.all.perms.h5
fi

PERM_FILES=$(echo fake.*.perms.h5)


#ANY PYTHON/2.7.* with h5py and numpy
module load enthought_python/7.3.2

#merge 10 perms at a time.
python ~/ESMtest/src/h5merge.py -i $PERM_FILES -o fake_merged.all.perms.h5 -n 10

