
#WTCCC Data Conversion for GWAS/ESM Test

WTCCC data came in two main types, the main genotype data and the supporting sample information file. In order to utilize these data in PLINK!, multiple preprocessing steps were taken. The goal being to convert the genotype/sample data into .ped/.map files which PLINK uses.
<br>

##Making the .map

The .map file is a necessary component for using PLINK!. The .map structure we are using is technically the "map3" format and it takes the following form:

|Chromsome #| rsID#|base pair position|
|:-----------|------------:|:-----------:|
|##|  rs# | ###### 

<br>


To make the .map, you will need the WTCCC genotype file, the Affymetrix 500k annotation, dbSNP SNP chromosome position table, and the dbSNP merger archive. Below, I will show how these are used and where to aqcuire them if neccessary.

<br>

###Gentoype files:
The genotype files are entitled, Affx_date_gt_disease_chiamo_chromosome.txt.gz. Each line in the genotype is an individuals genotype for a given SNP. Shown here is SNP id (rs or AFFY), individual ID (WTCCC), compound genotype(AA,GG,..etc), quality score.

zcat Affx_20070205fs1_gt_58C_Chiamo_22.txt.gz | head -n 5

|:-----------|------------:|:------------:|:------------:|
|rs915677|  WTCCC66751|	GG|	1.0000
|rs915677|	WTCCC66761|	GG|	1.0000
|rs915677|	WTCCC66771|	GG|	1.0000
|rs915677|	WTCCC66723|	GG|	1.0000
|rs915677|	WTCCC66731|	GG|	0.9971

<br>

However, notice that some SNPS are listed with Affymetrix ID:

zcat Affx_20070205fs1_gt_58C_Chiamo_22.txt.gz |cut -f 1 | sort | uniq | tail -n 5

|-----------:|------------:|
|SNP_A-4284676| | 
|SNP_A-4285289| |
|SNP_A-4293029| |
|SNP_A-4296471| |
|SNP_A-4298372| |  

<br>

###Affymetrix annotation files:

This is not directly useable because in order to use the most recent ENSEMBL dbSNP build to map SNPS to chromosome position we need all SNPS to be listed by rsID. So, the first step is to isolate and convert all SNP ID's to rsIDs. The necessary data for this can be found on the affymetrix product site http://www.affymetrix.com/support/technical/byproduct.affx?product=500k.The 500k SNP Chip is actually two sets of 250k Chips and thus the annotation files are separated by chip, namely you need: Mapping250K_Nsp Annotations, CSV format, Release 32 (83 MB, 7/15/11) and Mapping250K_Sty Annotations, CSV format, Release 32 (83 MB, 7/15/11). 

<br> 

To cut the list of SNPS from the genotype files, use the bash script getsnps.sh, which is centered around the command used to generatet the above table to create a file `chrom#.snps`. Then, send that list to the python program SNPID.py to create a fixed list of SNPS . This requires that you have downloaded the Affymetrix annotations, concatenated the Nsp and Sty files, and taken only the columns of Affy_id and Rs_ID into a separate .txt which SNPID.py will access. This can be done either programatically or with your favored spreadsheet application.

Make sure you change the line in SNPID.py where it references the affy annotation file:
```python
RS_file = io.open('/Path/to/annotation/AFFY_RSID','r', encoding='utf-16-le')
```

Then run is as such:
```bash
python /Path/to/SNPID.py 22.snps
```

<br>

###dbSNP tables:
Two dbSNP tables are necessary for making the .map: the SNP chromosome position table which is self explanatory and the merger archive which keeps track of SNP ID's which have been merged or changed over time. These can both be located on the dbSNP FTP server at ftp://ftp.ncbi.nih.gov/snp/database/organism_data/human_9606 and their specific names are b141_SNPChrPosOnRef.bcp.gz and RsMergeArch.bcp.gz . 

###Run the script to make the .map
Once you have the dbSNP tables and the processed list of snps, you can run the python script RSMerge.py which will take the dbSNP tables and list of snps and make the .map.

```sh
python /path/to/RSMerge.py /path/to/b141_SNPChrPosOnRef.bcp.gz /path/to/RsMergeArch.bcp.gz 22 DiseaseChrom#.map
```

This script will also output a file entitled Chrom#exclusion.txt, which will contain a list of rsID's which should be excluded from PLINK analysis on the basis of either missing data or not being in the Merge archive. While PLINK! had deal with missing data in normal .ped/.map format, it cannot do so in binary formats so we deal with missing data via exclusion list.

<br>

##Exclusion lists
In addition to the list of excluded SNP generated while making the .map, WTCCC has provided a list of SNPs for each disease which should be excluded for a variety of reasons. For each disease/chromosome combination the WTCCC and local exlusion lists should be concatenated to make a complete exclusion list. The WTCCC SNP exclusion lists need to processed to make sure all SNPs are properly identified by rsid. To do this, first you must run SNPID.py, as done with the total SNP list. Then you must run EXMerge.py, which is similiar to RSMerge.py in that it merges the SNP list with the newest dbSNP release.
<br.

##Making the .ped file

The .ped file is a necessary component for using PLINK!. The .ped structure we are using takes the following form:

|Fam | Indiv #|Pat|Mat| Sex|Pheno|Genotypes....|
|:-----------|------------:|:-----------|:-----------|------------:|:-----------:|:-----------:|
|0|  WTCCC# |0 |0|  0 |1/0|AA,GG,... 

We do not use any of the data except Indiv, Pheno and Genotypes, thus all other entries are set to 0. Technically, from PLINK's perspective we are using "compound genotypes"(AA vs A A), and 0/1 for phenotype corresponding to unaffected/affected.

<br>

To make the .ped, you will need the WTCCC genotype files, the WTCCC sample files, and the Affymetrix 500k annotation. The genotype and annotation files are explained above, they are not used in any different manner here. I will explain the use of the sample file.

###Sample Files:

The sample files contain information on the individuals in the panel. These sample files need to be screened agains the WTCCC Individual exclusion lists. A final filtered sample file for each disease should be used make the .ped file

###Run the script to make the .ped
Once you have the genotype and sample files located and the annotation file is properly made, you can use the python script WTCCCtoPed.py to create the .ped.

```sh
python /path/to/WTCCCtoPED.py /path/to/Affx_20070205fs1_gt_$Disease_Chiamo_$Chrom#.txt.gz $Chrom# /path/to/sample.txt /path/to/AFFY_RSID.txt $Disease /path/to/$Disease$Chrom#.ped
```

###Filtering and merging the .ped/.map for analysis
Before moving forward, it is best to use PLINK! to filter each .ped file for MAF<0.01. This will prevent the few SNPs which are messed up in just one disease panel from making it through to your final analysis. This can be done with the --recode option in PLINK:

```sh 

module load plink/1.90a

plink --file ${dis}${chrom}  --maf 0.01  --recode --out ${dis}${chrom}_filtered  
```

Then to merge your case and control map/ped files for an actual GWAS use the following command:

```sh 

module load plink/1.90a

plink --file ${dis}${chrom}_filtered --merge CONT${chrom}_filtered.ped CONT${chrom}_filtered.map --recode --out CONT_${dis}_${chrom} 
```


##Conversion to binary PLINK! format
Now that the .ped/.map files have been made we chose to use PLINK's inherent binary file capability for the sake of computation. PLINK itself can perform this operation when given a set of .ped/.map. This is a good time to filter the merged case control panel for genotyping rate and once again for MAF<0.01 

```sh

module load plink/1.90a

plink --file CONT_${dis}_${chrom} --exclude GW_exclude_snps_FINAL.txt --map3  --make-bed --out CONT_${dis}_${chrom}  --maf 0.01 --geno 0.01

```

This command will create a .bed/.bim file set which are the binary .ped/.map and a set of .fam/.nosex/.log which we will for the most part not be using.



