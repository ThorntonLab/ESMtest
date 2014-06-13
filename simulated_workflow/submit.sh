#!/bin/bash

#Submit to queues on UCI HPC that KRT lab has access on
qsub -q krt,bio,pub64 ~/src/ESMtest/simulated_workflow/00_makedata.sh $1
qsub -q krt,bio,pub64 ~/src/ESMtest/simulated_workflow/01_permute.sh $1
