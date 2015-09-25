# Small Example

The small example simulates 6000 individuals  with 500 SNPS. It is run
with one thread and can be simply executed with

cd small
sh submit.sh

This will make the following files
1. PLINK! format .ped/.bed ...etc files
2. PLINK! association test files
3. Permuted single marker test statistics in .h5 file format
4. ESM test p-values by loci midpoint (fake.esmpv.txt)

If you look in submit.sh, you will see each step. Let's go through
them

##makedata.sh:
	* ms simulates 1 locus with 500 SNPS for 6000 individuals with a
    recent severe bottleneck.
	*ms2plink converst ms output to PLINK! format
	*plink then converts to binary PLINK! format

##permute.sh:
	*Use a named pipe to handle permutation raw output from PLINK
	*plink does single marker test on permuted datasets. Here we do
	20,000 permutations in 2 sets
	*perms2h5 converts permutation output to h5 format for later use
	saving data in chunks of 50 markers by 10,000 perms

##esmk.sh:
	*Takes the h5 permutation files and runs the test
	*Window size = 10,000 BP,
	*Jump size = 1000 BP
	*Number of markers per window = 50
	*LD cutoff = 0.5
	*cmarkers 50 & cperms 10000 is the chunk size of the permutation
	data

Note that there are also merge.sh and esmk_merged.sh. Those scripts
are there incase you want to see how you would merge the h5 files and
run the test on the larger h5 file.
