#!/bin/bash

#Submit to KRT's queue on UCI HPC
qsub -q krt ~/src/ESMtest/simulated_workflow/00_makedata.sh $1
qsub -q krt ~/src/ESMtest/simulated_workflow/01_permute.sh $1
