dbSNP
====

We need to map rsID numbers to chromo and position.  The rdID numbers are from NCBI's dbSNP.

The numbers that we have for WTCCC are out of date. We can update them, though, with the latest dbSNP build.

You find dbSNP [here](http://www.ncbi.nlm.nih.gov/SNP/)

The ftp serve is [here](ftp://ftp.ncbi.nih.gov/snp/), and you can follow links to [here](ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/BED/), which is the human dbSNP info in BED format, which looks like this for chr1:

```{sh}
track name=dbSNP_human description="dbSNP Build 138 (GRCh37.p10)" date="2013-07-28 05:44" taxId=9606 dbSnpBuild=138 URL="http://www.ncbi.nlm.nih.gov/snp" assembly=GRCh37.p10 assemblyAccession=GCF_000001405.22
chr1    175261678       175261679       rs171   0       -
chr1    20869460        20869461        rs242   0       +
chr1    6160957 6160958 rs538   0       -
chr1    93617545        93617546        rs546   0       +
chr1    15546824        15546825        rs549   0       +
chr1    203713132       203713133       rs568   0       +
chr1    24181040        24181041        rs665   0       -
chr1    53679328        53679329        rs672   0       +
chr1    173876560       173876561       rs677   0       -
chr1    161191521       161191522       rs685   0       -
chr1    230845793       230845794       rs699   0       -
chr1    233971982       233971983       rs701   0       +
chr1    32372138        32372139        rs717   0       +
```

The first line really out to be commented out...

But, the columns are:

1. chromosome
2. start
3. stop
4. rsID number 
5. score.  often just left at 0
6. strand

Note: start and stop are with respect to position 1 being 0!!!!!  In the 1-offset world of molecular biology, that means that column 3 is the actual SNP position.
