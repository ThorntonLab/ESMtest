#!/bin/sh
#Make the data
sh makedata.sh
#Permute the data
sh permute.sh 
#Run the ESM test
sh esmk.sh


#OPTION TO MERGE THE H5FILES
#sh makedata.sh
#sh permute.sh
#sh merge.sh
#sh esmk_merged.sh



