#!/bin/bash

# Script to vuild trees using RAXML

cd /home/james/HT_Workflow/trees/Rex1
raxmlHPC -m GTRGAMMA -p 12345 -# 20 -s All_Rex1_filtered_mafft_gblocks.fasta -n T13 -o CR1-D
  
  raxmlHPC -m GTRGAMMA -p 12345 -b 12345 -# 100 -s All_Rex1_filtered_mafft_gblocks.fasta -n T14 -o CR1-D
  
  raxmlHPC -m GTRCAT -p 12345 -f b -t RAxML_bestTree.T13 -z RAxML_bootstrap.T14 -n T15 -o CR1-D


cd /home/james/HT_Workflow/trees/RTE
raxmlHPC -m GTRGAMMA -p 12345 -# 20 -s All_RTE_filtered_mafft_gblocks.fasta -n T13 -o L1-Snek_1
  
  raxmlHPC -m GTRGAMMA -p 12345 -b 12345 -# 100 -s All_RTE_filtered_mafft_gblocks.fasta -n T14 -o L1-Snek_1
  
  raxmlHPC -m GTRCAT -p 12345 -f b -t RAxML_bestTree.T13 -z RAxML_bootstrap.T14 -n T15 -o L1-Snek_1