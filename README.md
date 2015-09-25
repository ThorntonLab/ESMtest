#A PLINK-based implementation of the "ESM" test for associations due to rare alleles in GWAS data

Implements an association test based on the ESM statistic from [this paper](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1003258).

## Dependencies

1.  [GCC] (https://gcc.gnu.org/gcc-4.8/) -- version 4.8+
2.  [PLINK!](https://www.cog-genomics.org/plink2) -- 1.90a was used, but 1.07 also works.
3.  [boost](http://www.boost.org) --  Version 1.53 or greater is fine.
4.  [zlib](http://zlib.net) -- Version 1.2.7 is required.  mergeperms.cc checks this at compile time and will fail if a lower version number is encountered
5.  [GSL](http://gnu.org/software/gsl)
6.  [HDF5](https://www.hdfgroup.org/HDF5/release/obtain5.html) --
    version 1.8.11 or greater is fine (Install with --enable-cxx
    during configure step)
7.  [Python](https://www.python.org/downloads/)--2.7.2+ with [numpy](http://www.numpy.org/) and  [h5py](http://www.h5py.org/) -- Only needed if using h5merge.py,which is not strictly necessary


Please use your system's package installation tools to install the above whenever possible.

##Installing ESM


**git clone https://github.com/ThorntonLab/ESMtest.git ESMtest**

**cd ESMtest**

**./configure**

**make**

**make install**

If you want to install it in non-standard location (i.e your prefix if not /usr/local/) then:

**cd ESMtest**

**./configure --prefix=$HOME (or wherever you want it to go)**

**make**

**make install**
