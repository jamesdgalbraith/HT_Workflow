#!/bin/bash

cd ~/HT_Workflow/HT_Workflow/CARP/Aipysurus_laevis/long_seq

bundle -bundle 100000 -cut 600 -in long_carp_seqs.fasta

parallel -j 12 'FASTA={}; rpstblastn -evalue 0.01 -query $FASTA -db ~/Databases/localrpsb/transposons/transposons -num_threads 1 -out ${FASTA%.fa}.rps_out -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs stitle"' ::: *.fa

cat *out > long_carp_seqs_rps.tsv

rm long_carp_seqs.fasta-*
