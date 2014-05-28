#!/bin/bash

#Submit to KRT's queue on UCI HPC

for i in ~/src/ESMtest/simulated_workflow/0*.sh
do
    qsub -q krt $i ~/src/ESMtest
done