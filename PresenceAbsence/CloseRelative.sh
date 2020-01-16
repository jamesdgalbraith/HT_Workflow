#!/bin/bash

# Blast script used to identify highly similar LINEs in closely related elapids.
# Altered to search for each repeat in each species
# Example below  - searching for RTE-Snek in the brown snake genome
blastn -query RTE-Snek.fasta -db Genomes/Pseudonaja_textilis/EBS10Xv2-PRI.fasta -outfmt 6 | awk '{if ($3 > 97 && $4 > 100) print $0}' > RTE-Snek_in_EBS10Xv2-PRI_tsv