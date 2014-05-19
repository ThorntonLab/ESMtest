#!/bin/bash

PATH2PLINK=../../plink_test/plink-1.07-srcKRT

#Simulate some data.
ms 12000 1 -s 500 | ./ms2plink fake.ped fake.map
#Make binary input files
$PATH2PLINK/plink --noweb --file fake --make-bed --map3 --out fake
#Do perms
$PATH2PLINK/plink --noweb --bfile fake --assoc --map3 --mperm 1000000 --mperm-save-all --hwe 1e-6 --out fake &
gzip fake.mperm.dump.all