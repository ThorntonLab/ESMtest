#!/bin/sh

#Make the data
qsub  ~/ESMtest/simulated_workflow/00_makebigdata.sh

#Permute the data
qsub  ~/ESMtest/simulated_workflow/01_permutebig.sh

#OPTION 1 run the esm test on separate h5 files
qsub ~/ESMtest/simulated_workflow/02_esmkbig.sh

#OPTION 2 merge h5 files and run on big single file
#qsub ~/ESMtest/simulated_workflow/02_mergebig.sh
#qsub ~/ESMtest/simulated_workflow/03_esmk_mergedbig.sh 
