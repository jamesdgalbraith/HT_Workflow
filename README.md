# HT Workflow
Scripts used in the discovery, curation and analysis of horizonally transferred repeats in _Aipysurus laevis_

## Identification of repetitive sequences in _Aipysurus laevis_ (NovelLINEIdendification)
Repetitive sequences identified and classified using CARP (Zeng et al. 2018) and RepeatModeler (Smit and Hubley)

LINE_Finder.R - used to identify potential LINEs from Unclassified and Unknown repetitive sequences CARP and RepeatModeler based on protein domain presence

## Presence/absence of highly similar LINEs in closely related species
Command line BLASTN+ (megablast) other snake genomes for HTT LINEs. Also used to identify to HTT LINEs in _Aipysurus laevis_ transcriptome

## Presence/absence and curation of similar LINEs other highly divergent taxa

LINESearcher.R - used to identify species containing repeats similar to the HTT LINES

extendAlign.R and variants - used to construct alignments of similar sequences found in each species. Alignments manually edited and consensus made in Geneious

## Repeat phylogeny construction and figure building 
lines_for_trees.R - used to identify intact LINEs from RepBase and consensus sequences based on protein domain presence

mafft used to align intact LINEs

Trimming of alignment - online Gblocks with "Allow smaller final blocks", "Allow gap positions within the final blocks" and "Allow less strict flanking positions"

raxml.sh - used for tree building

## Repeat insertions near and in genes
overlapSearchConfirmations - used to identify insertions of HTT LINEs into near and into genes and confirm RepeatMasker annotation of said LINEs

insertionConfirmation.R - used to search for and create MSAs of ortholgous regions surrounding HTT LINEs inserted into and near genes in closely related species

