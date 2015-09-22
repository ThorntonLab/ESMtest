#!/bin/sh

if [ -a fake_merged.all.perms.h5 ]
then
    rm fake_merged.all.perms.h5
fi

PERM_FILES=$(echo fake.*.perms.h5)


#NEED:ANY PYTHON/2.7.* with h5py and numpy
module purge
module load krthornt/thorntonlab/1.0

#merge 10 perms at a time.
python ~/ESMtest/src/h5merge.py -i $PERM_FILES -o fake_merged.all.perms.h5 -n 10

