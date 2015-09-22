#A PLINK-based implementation of the "ESM" test for associations due to rare alleles in GWAS data

Implements an association test based on the ESM statistic from [this paper](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1003258).

## Dependencies

1.  [PLINK!](https://www.cog-genomics.org/plink2) -- 1.90a was used, but 1.07 also works.
2.  [boost](http://www.boost.org) --  Version 1.53 or greater is fine.
3.  [zlib](http://zlib.net) -- Version 1.2.7 is required.  mergeperms.cc checks this at compile time and will fail if a lower version number is encountered
4.  [GSL](http://gnu.org/software/gsl)
5.  [HDF5](https://www.hdfgroup.org/HDF5/release/obtain5.html) -- version 1.8.11 or greater is fine
6.  [h5py](http://www.h5py.org/) -- Only needed if using h5merge.py, which is not strictly necessary

Please use your system's package installation tools to install the above whenever possible.

## Dependencies for Simulated Examples

1. [ms](http://home.uchicago.edu/rhudson1/source/mksamples.html)
2. [libsequence](https://github.com/molpopgen/libsequence).

If you need to compile from source and are generally uncomfortable doing so, you may use this [script](https://github.com/molpopgen/install_libseq), which also installs [libsequence](https://github.com/molpopgen/libsequence).  Please read carefully the README that comes with the script.

##Installing ESM


**git clone https://github.com/ThorntonLab/ESMtest.git ESMtest**

**cd ESMtest**

**./configure**

**make**

**make install**

If you want to install it in non-standard location (i.e your prefix if not /usr/local/) then:

**cd ESMtest**

**./configure --prefix=$HOME (or wherever you want it to go**

**make**

**make install**
