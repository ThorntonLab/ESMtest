WTCCC data
=====

The data consist of 7 panels of cases plus two panels of controls.  The controls are merged into one big control group.

The data codes are:

1. 1958BC = 1958 Birth Cohort, control
2. NBS = National Blood Service Registry, control
3. BD = Bipolar disorder, case
4. CAD = Coronary Artery Disease, case
5. HT = Hypertension, case
6. IBD = Inflammatory Bowel Disease, case
7. RA = Rheumatoid Arthritis, case
8. T1D = Type I Diabetes, case
9. T2D = Type II Diabetes, case

The data files
====
The data files are in "[Chiamo](http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html)" format.  This format is intended for use with [SNPTEST](https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html), the software used for the original WTCCC analysis.  The output from SNPTEST is included with the data, but we've never looked at it.

On task that we need to accomplish is Chiamo-to-PLINK format conversion.  There are hints of tools available online:

1. [GenABEL](http://www.genabel.org/)

No idea if these work, or of the WTCCC formats are out-of-date with respect to what these files expect.

How are the data structured?

There are 25 gzipped text files with names like __Affx_20070205fs1_gt_HT_Chiamo_01.txt.gz__.  Files *\_01.txt.gz through *\_22.txt.gz represent the 22 autosomes.  The file ending in \_24.txt.gz appears to represent the X/Y psuedoautosomal region (this was inferred via Google searching the SNP IDs).  The files ending in \_23.txt.gz and \_X.txt.gz both appear to be X chromosome data (again, via Googling SNP IDs).  Curiously, the two files contain many of the same SNPs, but __NOT__ the same SNPs, which is totally bizarre:

```{sh}
zcat Affx_20070205fs1_gt_HT_Chiamo_23.txt.gz|cut -f 1 | sort | uniq > 23.snps
zcat Affx_20070205fs1_gt_HT_Chiamo_X.txt.gz|cut -f 1 | sort | uniq > X.snps
wc -l *.snps
 10274 23.snps
 10536 X.snps
 head -n 10 *.snps
==> 23.snps <==
rs1000489
rs1000530
rs1001874
rs1002116
rs1003752
rs10047021
rs1004823
rs1004991
rs1005155
rs1005303

==> X.snps <==
rs1000489
rs1000530
rs1001874
rs1002116
rs1003752
rs10047021
rs1004823
rs1004991
rs1005155
rs1005303
```

[This](http://www.wtccc.org.uk/info/data_formats.html) website makes the following statement regarding 23-vs-X:

"Note that, for some data sets on this site, the chromosome X data has been split into two 'chromosomes': 23 and 24. The region not homologous with Y (23) needed to be treated differently from the pseudo autosomal region (24)."

Ahh, standards...

Anyhow, we will probably follow tradition and ignore the X here anyways, focusing only in chr1-22.

__NOTE: the section below replaced real indentifiers, genotypes, etc. with REDACTED, as these data are covered by privacy concerns__

The genotype files have four columns.  These are SNP id, individual ID, genotype, and genotype score:

```{sh}
rs3094315  WTCCCREDACTED	REDACTED	1.0000
rs3094315	WTCCCREDACTED	REDACTED	1.0000
rs3094315	WTCCCREDACTED	REDACTED	1.0000
rs3094315	WTCCCREDACTED	REDACTED	1.0000
rs3094315	WTCCCREDACTED	REDACTED	1.0000
rs3094315	WTCCCREDACTED	REDACTED	1.0000
rs3094315	WTCCCREDACTED	REDACTED	1.0000
rs3094315	WTCCCREDACTED	REDACTED	1.0000
rs3094315	WTCCCREDACTED	REDACTED	0.9989
rs3094315	WTCCCREDACTED	REDACTED	1.0000
```

In addition, each panel has a "sample file" with lines like this:

```{sh}
WTCCCREDACTED      1       REDACTED      BRIGHT  11140A10        REDACTED  6       Unknown
WTCCCREDACTED      2       REDACTED      BRIGHT  11140A11        REDACTED        6       Unknown
WTCCCREDACTED      2       REDACTED      BRIGHT  11140A12        REDACTED        5       Unknown
```

The columns are:

1. Individual ID
2. Sex (male = 1, female = 2)
3. Panel label
4. Unknown
5. Unknown
6. Origin
7. Unknown
8. Unknown

Details of these file formats can be found [here](http://www.wtccc.org.uk/info/data_formats.html).




Exclusion lists
===
One complication of the data is that the raw genotype files cannot simply be used "as-is".  Each panel comes with "exclusion lists" for both individuals and SNPs:

```{sh}
find . -name "exclusion*"|grep -v gpg
./1958BC/20070205fs1/exclusion-list-05-02-2007.txt
./1958BC/20070205fs1/exclusion-list-snps-26_04_2007.txt
./1958BC/exclusion-list-05-02-2007.txt
./1958BC/exclusion-list-snps-26_04_2007.txt
./BD/exclusion-list-05-02-2007-BD.txt
./BD/exclusion-list-snps-26_04_2007.txt
./CAD/exclusion-list-05-02-2007-CAD.txt
./CAD/exclusion-list-snps-26_04_2007.txt
./HT/exclusion-list-snps-26_04_2007.txt
./HT/exclusion-list-05-02-2007-HT.txt
./IBD/exclusion-list-snps-26_04_2007.txt
./IBD/exclusion-list-05-02-2007-CD.txt
./NBS/20070205fs1/exclusion-list.txt
./NBS/20070205fs1/exclusion-list-snps-26_04_2007.txt
./NBS/exclusion-list.txt
./NBS/exclusion-list-snps-26_04_2007.txt
./RA/exclusion-list-snps-26_04_2007.txt
./RA/exclusion-list-05-02-2007-RA.txt
./T1D/exclusion-list-05-02-2007-T1D.txt
./T1D/exclusion-list-snps-26_04_2007.txt
./T2D/exclusion-list-05-02-2007-T2D.txt
./T2D/exclusion-list-snps-26_04_2007.txt
```

Let's look at a file excluding individuals.  We see four columns: panel label, SOMETHING, individual label reason for exclusion:

```{sh}
cat ./HT/exclusion-list-05-02-2007-HT.txt
#1. Affy CEL file bad
#2. Missing > 3%
#3. Het > 0.3 || Het < 0.225
#4a. External discordance: DIL genotyping
#4b. External discordance: NBS blood type
#4c. External discordance: phenotype measure in CD.
#5. Non-european ancestry
#6a. Sent twice by Sanger
#6b. Duplicated at unknown stage
#7. 1st/2nd degree relatives
REDACTED  19378A2	WTCCC169703	1,2
REDACTED	17947G12	WTCCC169139	2
REDACTED	18786F8	WTCCC169437	2
REDACTED	18786H5	WTCCC169489	2
REDACTED	18795B1	WTCCC169403	2
REDACTED	18795C6	WTCCC169367	2
REDACTED	18795G2	WTCCC169392	2
REDACTED	18796A7	WTCCC169977	2
#etc... (rest of file not shown)
```

And how the top of a SNP exclusion file.  The columns are pretty straightforward:

```{sh}
head -n 20 ./HT/exclusion-list-snps-26_04_2007.txt 
# The following 3 filters were applied sequentially. The number of the first filter a snp fails is gven in the FILTER column.
#1. (Studywise missing data proportion > 0.05) OR (Studywise minor allele frequency < 0.05 AND Studywise missing data proportion > 0.01)
#2. 58C+NBS HWE Exact Test p-value < 5.7e-7
#3. (58C vs NBS 1df Trend Test p-value < 5.7e-7) OR (58C vs NBS 2df General Test p-value < 5.7e-7)
CHR  AFFY_ID	RS_ID	FILTER
1	SNP_A-2169457	rs3817856	2
1	SNP_A-2218153	rs12141314	1
1	SNP_A-1842509	rs2292857	1
1	SNP_A-4218776	rs262683	2
1	SNP_A-2038545	rs4648515	1
1	SNP_A-1785968	rs12119163	2
1	SNP_A-4242198	rs4648451	1
1	SNP_A-4250429	rs2981884	1
1	SNP_A-2242875	rs6704012	1
1	SNP_A-2107218	rs2651912	1
1	SNP_A-4218510	rs6670518	1
1	SNP_A-2206691	rs6690558	1
1	SNP_A-2219183	rs2483280	1
1	SNP_A-1786405	rs16823663	1
1	SNP_A-2258506	rs13375075	1
```

So what needs to be done?
====

We need the following:

1.  Chiamo-to-PED/MAP conversion.  It looks like there are tools to take PED/MAP to other formats, so we may have to write our own.
2.  It is almost certainly the case that the MAP files must be made separately for each case/control pair.  __We should reread the WTCCC paper and check if it states whether or not markers/individuals were excluded globally or on a per-comparison basis.  I suspect that it is the latter.__
3.  The exclusion files need to be accounted for while converting the Chiamo files.  Exactly how depends on how we treat the issues brought up in the previous point.

This should result in 7 pairs of PED/MAP files, residing in subdirectories named using the data codes described above.  We can then use PLINK for remaining filtering (HWE, LD, etc.)

One important point:  the exclusion lists describe markers that are excluded due to HWE violations.  However, my own previous analysis found that there are more markers violating HWE.  We'll use PLINK to handle that.

It looks like we really should use PLINK's utility to convert from PED to binary (see PLINKnotes).  At the conversion step, we can exclude individuals and snps.  It also looks like we can filter on HWE, etc., as well. The HHE filtering during --make-bed is not documented that I can see, but I tested it on the example data:

```{sh}
../plink-1.07-src/plink --file ../example/wgas1 --make-bed --hwe 0.05

@----------------------------------------------------------@
|        PLINK!       |     v1.07      |   10/Aug/2009     |
|----------------------------------------------------------|
|  (C) 2009 Shaun Purcell, GNU General Public License, v2  |
|----------------------------------------------------------|
|  For documentation, citation & bug-report instructions:  |
|        http://pngu.mgh.harvard.edu/purcell/plink/        |
@----------------------------------------------------------@

Web-based version check ( --noweb to skip )
Recent cached web-check found... OK, v1.07 is current

+++ PLINK 1.9 is now available! See above website for details +++ 

Writing this text to log file [ plink.log ]
Analysis started: Fri May 16 15:47:01 2014

Options in effect:
  --file ../example/wgas1
	--make-bed
	--hwe 0.05

228694 (of 228694) markers to be included from [ ../example/wgas1.map ]
90 individuals read from [ ../example/wgas1.ped ] 
90 individuals with nonmissing phenotypes
Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)
Missing phenotype value is also -9
49 cases, 41 controls and 0 missing
45 males, 45 females, and 0 of unspecified sex
Before frequency and genotyping pruning, there are 228694 SNPs
90 founders and 0 non-founders found
4579 markers to be excluded based on HWE test ( p <= 0.05 )
	5225 markers failed HWE test in cases
	4579 markers failed HWE test in controls
Total genotyping rate in remaining individuals is 0.993346
0 SNPs failed missingness test ( GENO > 1 )
0 SNPs failed frequency test ( MAF < 0 )
After frequency and genotyping pruning, there are 224115 SNPs
After filtering, 49 cases, 41 controls and 0 missing
After filtering, 45 males, 45 females, and 0 of unspecified sex
Writing pedigree information to [ plink.fam ] 
Writing map (extended format) information to [ plink.bim ] 
Writing genotype bitfile to [ plink.bed ] 
Using (default) SNP-major mode

Analysis finished: Fri May 16 15:47:11 2014

```
