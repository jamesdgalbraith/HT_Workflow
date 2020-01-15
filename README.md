# HT Workflow
Scripts used in the discovery, curation and analysis of horizonally transferred repeats in _Aipysurus laevis_

## Identification of repetitive sequences in _Aipysurus laevis_
Repetitive sequences identified and classified CARP (Zeng et al. 2018) and RepeatModeler (Smit and Hubley)
LINEs identified from Unclassified and Unknown repetitive sequences CARP and RepeatModeler using LINE_Finder.R

## Presence/absence of highly similar LINEs in closely related species
BLAST other snake genomes for LINEs/BLAST snake transcriptomes for LINEs - bash

## Presence/absence and curation of similar LINEs other highly divergent taxa
repeat_blaster.R
extendAlign.R
extendAlignSolo.R

### Tree building
Intact LINEs itendified with LINE_finder.R
Alignment of intact LINEs - bash mafft
Trimming - online Gblocks with "Allow smaller final blocks", "Allow gap positions within the final blocks" and "Allow less strict flanking positions"
Tree building - RaXML