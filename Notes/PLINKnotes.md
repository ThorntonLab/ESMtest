Using PLINK to permute GWAS data
======

[PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/) is a platform for analysing GWAS data.  Lots of features, etc.

Ineresting fact.  There is a "[fork](https://www.cog-genomics.org/plink2)" of PLINK that is a collab b/w the original authors, BGI, and some other folks.  It appears to be geared more towards nextgen-based GWAS. We'll stick with 1.07 for now...

One possible use of PLINK would be to implement the ESM test using PLINK to do all the heavy lifting.  This document is a set of notes on the good and bad of PLINK, with special emphasis on permutation output vis-a-vis best practices on high-performance computing clusters.

Compilation Notes
===
I have installed v1.07 from source on my iMac (OS X Mavericks, clang compilers).  It needed modification to several source files to get rid of compiler errors.  Also, the "SYS" variable needs to be set in the Makefile.  Instructions are provided in the Makefile.  Default is for Unix, and I had to set it to "MAC".

Test data
===
I downloaded "example.zip" from [here](http://pngu.mgh.harvard.edu/~purcell/plink/res.shtml#teach).  It contains the following stuff:

```{sh}
ls -lhrt example/
total 196120
-rw-r--r--@ 1 krthornt  staff   6.7M May 24  2008 wgas1.map
-rw-r--r--@ 1 krthornt  staff   1.6K May 24  2008 pop.cov
-rw-r--r--@ 1 krthornt  staff   1.7M May 24  2008 gPLINK.jar
-rw-r--r--@ 1 krthornt  staff   509B May 24  2008 extra.map
-rw-r--r--@ 1 krthornt  staff   5.0M May 24  2008 Haploview.jar
-rw-r--r--@ 1 krthornt  staff    79M May 24  2008 wgas1.ped
-rw-r--r--@ 1 krthornt  staff   8.1K May 26  2008 extra.ped
-rw-r--r--@ 1 krthornt  staff   3.8M May 26  2008 plink.exe
-rw-r--r--@ 1 krthornt  staff   2.7K Jul 11  2008 command-list.txt
```

The file "wgas1.ped" is the main input file.  The PED format is described [here](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped).

Each PED file must have a corresponding MAP file, which is info about each SNP.  Details are [here](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#map).

PLINK formatting issues
===
There are several issues for how to write the data files that we need to address:


1. For case/control data, 1 = control/unaffected, 2 = case/affected, 0 = missing.
2. We probably do not want to have sex identifiers, etc., in the PED files.
3. MAP files allow you to put in a negative number for the mutation position, indicating an excluded SNP.  This suggests the possibility of pre-filtering (based on HWE p-values, minimum allele counts, etc.), making a new MAP file representing the filtering, and then permuting.  However, the documentation says that this only works for plain-text input, not the binary-format input (which is desirable as it can be read into memory faster).  That may not matter--maybe the plain text is ok?

PLINK usage questions/comments
====

Questions:

1. Does PLINK have a way to subset the data?  For example, can we take large file and break it into chunks of $M$ markers?  it [this](http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#tabset) relevant?  Upon further poking around, [this](http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#extract) is probably what we are looking for.  PLINK allows us to get all SNPs in a window, meaning we can break the genome up into windows of a constant size, and permute the SNPs w/in that window.

Comments:

1. It [looks like](http://pngu.mgh.harvard.edu/~purcell/plink/thresh.shtml#hwd) we can exclude SNPs based on HWE in controls "on the fly", which is a huge simplification.
2. It is critical to keep the following in mind: we wish to permute SNPs in windows sizes that are chosen for computational convenience.  However, we will slide windows of __a different size__ over the data in order to calculate ESM.  This means that the permutations __must be held constant, genome-wide__, a la [Churchill and Doerge](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1206241/pdf/ge1383963.pdf).  If we fail to do this, then windows used in the ESM calculations will be mixtures of different underling permutation samples.  In practice, this means that all perms from the same chromosome must use the __same random number seed__.

Basic usage
=====

Maybe a useful [reference](http://128.40.230.201/stupkalab/downloads/2010.02.17_Lescai_TutorialPLINK.pdf)?

Vanilla association test.  The output is a $\chi^2$ statistic for each marker:
```{sh}
#!sh

#Basic association test, HWE test, and MAF counting, all in 1 command.  That is pretty cool.
#The --counts option replaces allele frequencies with counts in the output
../plink-1.07-src/plink --file ../example/wgas1 --assoc --counts --hardy --freq
```

Generating permutations of the data
===
PLINK can generate permutations of the data, and output either "empirical" p-values or the permuted values of all test statistics.

Several examples:
```{sh}
#!sh

#Exploring permutations.

#Generate 100 perms of the example data, with association statistics written out
#This command simply reports two types of permutation p-value for the association statistic
#../plink-1.07-src/plink --file ../example/wgas1 --assoc --mperm 100 --seed 101

#Let's save ALL the perms.
#The output from this is UGLY.  For K snps:
#perm number permstat1 permstat2 ... permstat_k
#The command below is the preferred input for ESM, as the chi-squared is rapid to calculate.  
#Although the p-value of ANY single-marker test can be used as the ESM input, 
#other tests can be 10-100x slower to compute PER MARKER, which is really, really bad when you are permuting
../plink-1.07-src/plink --file ../example/wgas1 --assoc --mperm 10000 --seed 101 --mperm-save-all

#We can do it w/logistic test.  This is SLOW, as expected
../plink-1.07-src/plink --file ../example/wgas1 --logistic --mperm 100 --seed 101 --mperm-save-all

#We can do it w/Fisher's Exact test.  Again, slow
../plink-1.07-src/plink --file ../example/wgas1 --fisher --mperm 100 --seed 101 --mperm-save-all
```

On my iMac, it took 37 minutes and 50 seconds to do the $10^4$ permutations of the example data, which consists of 90 individuals and over $2 \times 10^5$ SNPS:

```{sh}
wc example/wgas1.ped 
      90 41165460 82332090 example/wgas1.ped
wc testing/plink.frq.count 
  228695 1600865 10977360 testing/plink.frq.count
```

That timing isn't so bad, but it would be useful to test how long it takes to do, say, $3 \times 10^6$ permutations of $3\times 10^3$ individuals for 50 or 100 SNPs.  Presumably PLINK is smart enough that the $\chi^2$ calculation is the slow part.  In other words, I'm hoping that it is permuting indexes representing each individual rather than actually moving all the genotype data around.

The output file is 17GB, gzipping down to 5.3GB, and the gzipping took a LONG time.

Some issues immediately become apparent:

1. PLINK writes output to the screen.  In a GE environment, we'd need to redirect stdout and stderr to /dev/null, otherwise we'd crush the file system with the .e/.o files.
2. PLINK is not buffering output internally.  Rather, it continuously writes to a file called "plink.mperm.dump.all", which is plain-text.  This approach has several drawbacks.  First, it is not clear that a file system like gluster could handle this.  
3. PLINK has a set of default output file [names](http://pngu.mgh.harvard.edu/~purcell/plink/reference.shtml#output).  The default behavior is to name each output file plink.suffix, where suffix is a function of the arguments you are using.  The prefix may be changes via the --out option (see example somewhere down below).

The reason why these issues are important is that we will need to split up the data into chunks of $M$ markers, and generate a large number of permutations per chunks.  Because we cannot control the output name, we may need to make a subdirectory for each chunk.  Further, the uncompressed files are a problem in terms of space, but we can easily gzip after each run.

Some solutions for dealing w/the output include:

1. gzipping immediately after the perms are generated.  This is a no-brainer, but time-intensive.
2. adding the output files to a growing tar-format archive.  This requires file-locking at the shell level, which I have no direct experience with.  In principle, these tar archives are "mountable" as if they were their own partition.  Again, no experience, and may have to ask Harry.
3. Can we use mkfifo to simultaneously permute and gzip?

Yes, the named pipes do seem to work, which is great:
```{sh}
#!sh

#Exploring permutations using named pipes

mkfifo plink.mperm.dump.all
../plink-1.07-src/plink --file ../example/wgas1 --assoc --mperm 100 --seed 101 --mperm-save-all &
cat plink.mperm.dump.all | gzip > plink.mperm.dump.all.gz
rm -f plink.mperm.dump.all
```

No idea what the effect of the named pipe technique would be on gluster load, etc.

Anyhow, let's redo the $10^4$ perms of the test data using this technique:

```{sh}
#!sh

mkfifo plink.mperm.dump.all
../plink-1.07-src/plink --file ../example/wgas1 --assoc --mperm 10000 --seed 101 --mperm-save-all &
cat plink.mperm.dump.all | gzip > plink.mperm.dump.all.gz
rm -f plink.mperm.dump.all
```

Executing the script gave the following results:
```{sh}
time sh permutation2.sh

real  46m10.570s
user	19m55.896s
sys	0m38.645s
```

So, that's a 20% performance hit, but the benefit is that we skip the generation of one very large file.

Permuting subsets of the data
====
Let's try the following:

1.  Permute just a subset of the test data
2.  Apply a HWE filter.  The filter will be liberal, so that SNPs are certain to get tossed.
3.  Manipulate the output file names

```{sh}
#!sh

OFILEBASE=subset
OFILE=$OFILEBASE.mperm.dump.all
mkfifo $OFILE
#The arguments are nicely self-explanatory...
../plink-1.07-src/plink --file ../example/wgas1 --assoc --mperm 10000 --seed 101 --mperm-save-all --chr 2 --from-kb 5000 --to-kb 10000 --hwe 0.05 --out $OFILEBASE &
cat $OFILE | gzip > $OFILE.gz
rm -f $OFILE
```

Some of the info printed to stdout is encouraging, suggesting that SNPs are excluded based on rejected HWE in controls only (which is what the documentation says, so yay for that):
```{sh}
5 markers to be excluded based on HWE test ( p <= 0.05 )
  11 markers failed HWE test in cases
	5 markers failed HWE test in controls
  ```
  
  Additional notes:
  
  1. If you input --to-kb values that are greater than the length of the chromosome, the program handles it fine, and simply does the requested operation up until the last SNP.
  2. If your --from-kb is greater than the position of the last marker, PLINK just prints an error message to stderr

  Trying to tie it all together
  ====
  
We have figured out how to:
  
  1.  Calculate $\chi^2$ for every marker
  2.  Subset the data based on physical ranges
  3.  Obtain the permutation distribution of the $\chi^2$ statistic for the subsetted data
  4.  Use named pipes to avoid writing the very large permutation files.
  5.  How to control the prefix of the output file names
  6.  If you pick a window where there are no markers, the program exist gracefully.
  
When the run is successful, the output looks like:
  ```{sh}
   ls -lhrt
total 128
-rw-r--r--+ 1 krthornt  staff   324B May 15 13:08 permutation2.sh
-rw-r--r--+ 1 krthornt  staff    46K May 15 13:09 subset.mperm.dump.all.gz
-rw-r--r--+ 1 krthornt  staff   2.0K May 15 13:09 subset.log
-rw-r--r--+ 1 krthornt  staff   264B May 15 13:09 subset.assoc.mperm
-rw-r--r--+ 1 krthornt  staff   576B May 15 13:09 subset.assoc
  ```

The files are:

1. subset.log: a repeat of what is written to stdout.  This can be deleted, or they can all be archived for posterity's sake.
2. subset.assoc.mperm: the permutation p-value of the association test statistics.  We do not need this--delete
3. subset.assoc: the output of the association test.  This file is needed, as it contains the names of the SNPs that were analyzed, and is critical info, allowing us to read the MAP file back in and assign chromosome, position, etc.
4. subset.mperm.dump.all.gz: for each permuted data set, you get the permutation ID number, and then the permuted association test statistic for each marker, in the order that they appear in subset.assoc.  (That last statement is an assumption--not verified!)

I think that the situation for actually doing the perms is simple:

1.  Choose a window size corresponding to a number of SNPs that is modest on average (50? 100?).
2.  The output file name prefix can simply be chrC.start.stop, where C is the chromosome label (1-22, as the data are only autosomal), and the start/stop are the values passed to PLINK's --from-kb and --to-kb arguments.

The appeal is that we need very few tools (see next section).  The caveat is whether or not HPC can handle the I/O load.  That's Harry's problem...

So what do we need?
====

1.  A program to convert WTCCC data into PLINK format
2.  A program to merge the permutation output from PLINK back into a large file.  Given the sizes, I think we need to suck it up and write directly to gzipped files using the low-level [zlib](http://zlib.net/) library, which allows appending, seeking, etc.  This file needs to be "sensible", containing the SNP names from the WTCCC input files, so that they can be mapped directly to genes, etc., using current annotations.  Not doing this the first time through was a big impediment to downstream analysis.   In practice, this step merges the .assoc file with the .mperm.dump.all.gz file with the MAP file for the appropriate chromosome.
3.  A program to perform the ESM test on the output of #2
4.  We need to deal with LD between SNPS.  Not sure if this is a custom program to write, or if PLINK can do something automatically.

Outside-the-box ideas:

1.  Instead of piping the named pipe to gzip, could we simply buffer the data internally and then write it to a tar file?  Problem: large data = risk of node crashing.

Dealing with linkage disequilibrium (LD)
====

PLINK has [functions](http://pngu.mgh.harvard.edu/~purcell/plink/ld.shtml) for calculating LD.  They are easy to use.

To calculate pairwise correlation coefficients ( $r^2$ ) for __ALL__ pairs of markers:

```{sh}
../plink-1.07-src/plink --file ../example/wgas1 --r2
```

For just a single chromosome, outputting to chr2.ld and removing markers with significant HWE deviations in controls:

```{sh}
../plink-1.07-src/plink --file ../example/wgas1 --r2 --chr 2 --out chr2 --hwe 0.05
```

Sneakily, PLNK has some defaults that we need to over-write.  The following command calculates LD for all pairs of SNPs within 1000 SNPs of each other, within the same 10 megabase window (10,000kb), and only reports LD statistics > 0.8:

```{sh}
../plink-1.07-src/plink --file ../example/wgas1 --r2 0.8 --chr 2 --out chr2 --hwe 0.05 --ld-window 1000 --ld-window-kb 10000 --ld-window-r2 0.8
```

So, we just need to run that per chromosome once, gzip up the output, and store it.

Random notes about PLINK
====

Various insights copy/pasted from their website:

1. If an input file is compressed (gzip compression) and ends in the .gz extension, PLINK will automatically decompress it (if compiled with ZLIB support).  This appears to only apply to the [meta analysis module](http://pngu.mgh.harvard.edu/~purcell/plink/metaanal.shtml).  Lame.



