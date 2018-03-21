#!/usr/bin/bash

#SBATCH --nodes 1 --ntasks 8 --mem 8G --time 8:00:00 -p batch
CPU=8
module unload python
module load qiime
core_diversity_analyses.py -i Sam1_34a.otu_table.taxonomy.biom -m 050517NP515F-mapping2_NP.txt -o core_diversity -a --sampling_depth 5000 --tree_fp Sam1_34a.tree.phy  --recover_from_failure -O $CPU -c VegZone,Disturbance



