#Examples run on simulated data

Implments the whole ESM test workflow on data simulated in ms.

You must be aware of where/how the dependencies were installed. Are
you dependencies installed in the "usual places"?

## Small Example
500 SNPs for 6000 individuals

###permute data:
       *plink does single marker test on permuted datasets. Here we do
        2,000 permutations in 2 sets
plink --noweb --file fake --map3 --r2 --out fake.1 

plink --noweb --file fake --assoc --map3 --mperm 1000 --mperm-save-all  --out fake.1 --seed 1 

plink --noweb --file fake --map3 --r2 --out fake.2 

plink --noweb --file fake --assoc --map3 --mperm 1000 --mperm-save-all  --out fake.2 --seed 1 

      *perms2h5 converts permutation output to h5 format for later use
        saving data in chunks of 50 markers by 10,000 perms

perms2h5 -i fake.1.mperm.dump.all -o fake.1.perms.h5 -b fake.bim -n 50 -l fake.1.ld

perms2h5 -i fake.2.mperm.dump.all -o fake.2.perms.h5 -b fake.bim -n 50 -l fake.1.ld

rm -f fake.*.mperm.dump.all

rm -f fake.*.assoc.mperm

###run the esm test:
	*Takes the h5 permutation files and runs the test
	*Window size = 10,000 BP,
	*Jump size = 1000 BP
	*Number of markers per window = 50
	*LD cutoff = 0.5
	*cmarkers 50 & cperms 1000 is the chunk size of the permutation
	data

esmk -o fake.esmpv.txt -w 10000 -j 1000 -k 50 -n 1 -r 0.5 --cmarkers 50 --cperms 1000 --nperms 2000 fake.1.perms.h5 fake.2.perms.h5

Note that there are also merge.sh and esmk_merged.sh. Those scripts
are there incase you want to see how you would merge the h5 files and
run the test on the larger h5 file.


