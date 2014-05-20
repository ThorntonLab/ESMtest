#!/bin/bash

#Problem: this brakes the mkfifo procedure!!

PATH2PLINK=../../plink_test/plink-1.07-srcKRT

#Simulate some data.
ms 12000 1 -s 500 | ./ms2plink testdata.ped testdata.map
#Make binary input files
$PATH2PLINK/plink --noweb --file testdata --make-bed --map3 --out testdata
#Do perms
$PATH2PLINK/plink --noweb --bfile testdata --assoc --map3 --mperm 100000 --mperm-save-all --hwe 1e-6 --out testdata 

