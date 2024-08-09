#!/bin/bash

# #$ -M wstamey@nd.edu
# #$ -m abe
#$ -q long
#$ -N meta2x2t1n2
#$ -t 1-16

n_iter="1000"

cd "/afs/crc.nd.edu/user/w/wstamey/Private/Research Projects/meta2x2/tablesandgraphs/table1and2job/"
module load jags
module load R
Rscript ../paper_table_1and2_sims_trinom.R "$n_iter" "$SGE_TASK_ID"
