#!/bin/bash

#$ -M wstamey@nd.edu
#$ -m abe
#$ -q long
#$ -N meta2x2t1n2
module load jags
module load R
Rscript paper_table_1and2_sims.R
