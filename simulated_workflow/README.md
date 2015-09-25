#Examples run on simulated data

Implments the whole ESM test workflow on data simulated in ms.

You must be aware of where/how the dependencies were installed. Are
you dependencies installed in the "usual places"?

## Small Example
500 SNPs for 6000 individuals

###permute data:

*Use a named pipe to handle permutation raw output from PLINK
mkfifo fake.1.mperm.dump.all
mkfifo fake.2.mperm.dump.all

*plink does single marker test on permuted datasets. Here we do
        20,000 permutations in 2 sets

plink --noweb --bfile fake --r2 --assoc mperm=10000 --mperm-save-all --map3 --out fake.1 --seed 1 --threads 1 &
plink --noweb --bfile fake --r2 --assoc mperm=10000 --mperm-save-all --map3 --out fake.2 --seed 2 --threads 1 &

*perms2h5 converts permutation output to h5 format for later use
        saving data in chunks of 50 markers by 10,000 perms

perms2h5 -i fake.1.mperm.dump.all -o fake.1.perms.h5 -b fake.bim -n 50 -l fake.1.ld
perms2h5 -i fake.2.mperm.dump.all -o fake.2.perms.h5 -b fake.bim -n 50 -l fake.1.ld
rm -f fake.*.mperm.dump.all
rm -f fake.*.assoc.mperm

###esmk.sh:
	*Takes the h5 permutation files and runs the test
	*Window size = 10,000 BP,
	*Jump size = 1000 BP
	*Number of markers per window = 50
	*LD cutoff = 0.5
	*cmarkers 50 & cperms 10000 is the chunk size of the permutation
	data
esmk -o fake.esmpv.txt -w 10000 -j 1000 -k 50 -n 1 -r 0.5 --cmarkers 50 --cperms 10000 --nperms 20000 fake.1.perms.h5 fake.2.perms.h5

Note that there are also merge.sh and esmk_merged.sh. Those scripts
are there incase you want to see how you would merge the h5 files and
run the test on the larger h5 file.


## Big Example for Simulated Examples

The big example is the same in many ways the small example except that
it simulated 40000 SNPs and perfomrs 100K permutations. The scripts
are thus structured for submission to cluster queuing system. This is
not necessary and you can modify them to work on any system with 4 or
more cores.

