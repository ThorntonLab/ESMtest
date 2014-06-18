#!sh
#$ -q krt,bio,pub64
#$ -N GETSNPS

for i in /bio/krthornt/WTCCC/WTCCC_Jan2013/1958BC/20070205fs1/*Chiamo*.gz
do
a=$(echo $i |rev | cut -d_ -f1 | rev);
b=$(echo $a |rev | cut -d. -f3 | rev);

zcat $i |cut -f 1 | sort | uniq > $b.snps
done 
 
