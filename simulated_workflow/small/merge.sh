#!/bin/sh

if [ -a fake_merged.all.perms.h5 ]
then
    rm fake_merged.all.perms.h5
fi

PERM_FILES=$(echo fake.*.perms.h5)


#On UCI HPC we use                
#module purge
#module load krthornt/thorntonlab/1.0
#WHICH LOADS:
##gcc/4.8.4
##gsl/1.16
##htslib/1.2.1
##zlib/1.2.8
##boost/1.59.0
##hdf5/1.8.15-patch1
##python/2.7.10
##krthornt/libsequence/1.8.5

#module load plink/1.90a

#Modify as needed for your system                                                                                                                 
#NEED:
#python 2.7+ with h5py and numpy

#merge 10 perms at a time.
python ~/ESMtest/src/h5merge.py -i $PERM_FILES -o fake_merged.all.perms.h5 -n 10

