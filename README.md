# HT Workflow
Scripts used in the discovery, curation and analysis of horizonally transferred repeats in _Aipysurus laevis_

## Identification of repetitive sequences in _Aipysurus laevis_
Repetitive sequences identified and classified using CARP (Zeng et al. 2018) and RepeatModeler (Smit and Hubley)

LINE_Finder.R - used to identify potential LINEs from Unclassified and Unknown repetitive sequences CARP and RepeatModeler based on protein domain presence

## Presence/absence of highly similar LINEs in closely related species
BLAST other snake genomes for LINEs/BLAST snake transcriptomes for LINEs

## Presence/absence and curation of similar LINEs other highly divergent taxa
repeat_blaster.R - used to identify similar LINEs using relaxed BLAST searchs

extendAlign.R and variants - used to construct alignments of similar sequences found in each species. Alignments manually edited and consensus made in Geneious

## Figure building - tree building and consensus plotting
lines_for_trees.R - used to identify intact LINEs from RepBase and consensus sequences based on protein domain presence

mafft used to align intact LINEs

Trimming of alignment - online Gblocks with "Allow smaller final blocks", "Allow gap positions within the final blocks" and "Allow less strict flanking positions"

raxml.sh - used for tree building