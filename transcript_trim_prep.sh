#!/bin/bash

# script used to prepare HT repeats for trimming of 3' end

# search for transcript, if above 99% id and 1000 base pair hit print into bed format
blastn -query Proto2-Snek_transcript.fa -db ~/Analysis/Genomes/Aipysurus_laevis/kmer_49.pilon_x2.sorted.fasta -outfmt 6 | awk '{if ($3 > 99 && $4 > 1000 && $10 > $9) print $2 "\t" $9-1 "\t" $10 "\t.\t.\t+"; if ($3 > 99 && $4 > 1000 && $10 < $9) print $2 "\t" $10-1 "\t" $9 "\t.\t.\t-"}' > Proto2-Snek_transcript_hits.bed

# get fasta of repeat hits
bedtools getfasta -s -fi ~/Analysis/Genomes/Aipysurus_laevis/kmer_49.pilon_x2.sorted.fasta -fo Proto2-Snek_transcript_hits.fa -bed Proto2-Snek_transcript_hits.bed 

# combine initial query and transcript hits
cat Proto2-Snek_transcript.fa Proto2-Snek_transcript_hits.fa > Proto2-Snek_transcript_hits_compiled.fa

# align using mafft
mafft --localpair Proto2-Snek_transcript_hits.fa > Proto2-Snek_transcript_aligned_hits.fa

# manually trim transcript in Geneious
# used manually trimmed transcripts in place of consensus going forward