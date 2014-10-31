#Tools for processing PLINK 1.90a's permutation output

##Dependencies

1.  [boost](http://www.boost.org) --  Version 1.53 or greater is fine.
2.  [zlib](http://zlib.net) -- Version 1.2.7 is required.  mergeperms.cc checks this at compile time and will fail if a lower version number is encountered
3.  [GSL](http://gnu.org/software/gsl)
4.  [h5py](http://www.h5py.org/) -- Only needed if using h5merge.py, which is not strictly necessary

Please use your system's package installation tools to install the above whenever possible.

If you need to compile from source and are generally uncomfortable doing so, you may use this [script](https://github.com/molpopgen/install_libseq), which also installs [libsequence](https://github.com/molpopgen/libsequence).  Please read carefully the README that comes with the script.

##Assumptions