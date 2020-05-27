library(tidyverse)
library(Biostrings)

# set query and genome
query_repeat <- "~/ElapidRepeats/CARP/Aipysurus_laevis/long_seq/potential_LTR_elements.fa"

# perform blastn search
blastn_aipysurus <- read_tsv(system(paste0("blastn -num_threads 12 -dust yes -query ", query_repeat, " -db ~/Genomes/Reptiles/Aipysurus_laevis/kmer_49.pilon_x2.sorted.fasta -outfmt \"6 qseqid sseqid qstart qend sstart send pident qcovs bitscore length mismatch evalue qlen slen\""), intern = TRUE), col_names = c("qseqid", "seqnames", "qstart", "qend", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "qlen", "slen")) %>%
  as_tibble() %>%
  dplyr::filter(length >= 0.2 * qlen)

aipysurus_table <- tibble(qseqid = names(table(blastn_aipysurus$qseqid)), aipysurus_n = table(blastn_aipysurus$qseqid))

print("Aipysurus_laevis")

blastn_emydocephalus <- read_tsv(system(paste0("blastn -num_threads 12 -dust yes -query ", query_repeat, " -db ~/Genomes/Reptiles/Emydocephalus_ijimae/emyIji_1.0.fasta -outfmt \"6 qseqid sseqid qstart qend sstart send pident qcovs bitscore length mismatch evalue qlen slen\""), intern = TRUE), col_names = c("qseqid", "seqnames", "qstart", "qend", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "qlen", "slen")) %>%
  as_tibble() %>%
  dplyr::filter(length >= 0.2 * qlen)

emydocephalus_table <- tibble(qseqid = names(table(blastn_emydocephalus$qseqid)), emydocephalus_n = table(blastn_emydocephalus$qseqid))

print("Emydocephalus_ijimae")

blastn_hydrophis <- read_tsv(system(paste0("blastn -num_threads 12 -dust yes -query ", query_repeat, " -db ~/Genomes/Reptiles/Hydrophis_melanocephalus/hydMel_1.0.fasta -outfmt \"6 qseqid sseqid qstart qend sstart send pident qcovs bitscore length mismatch evalue qlen slen\""), intern = TRUE), col_names = c("qseqid", "seqnames", "qstart", "qend", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "qlen", "slen")) %>%
  as_tibble() %>%
  dplyr::filter(length >= 0.2 * qlen)

hydrophis_table <- tibble(qseqid = names(table(blastn_hydrophis$qseqid)), hydrophis_n = table(blastn_hydrophis$qseqid))

print("Hydrophis_melanocephalus")

blastn_notechis <- read_tsv(system(paste0("blastn -num_threads 12 -dust yes -query ", query_repeat, " -db ~/Genomes/Reptiles/Notechis_scutatus/TS10Xv2-PRI.fasta -outfmt \"6 qseqid sseqid qstart qend sstart send pident qcovs bitscore length mismatch evalue qlen slen\""), intern = TRUE), col_names = c("qseqid", "seqnames", "qstart", "qend", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "qlen", "slen")) %>%
  as_tibble() %>%
  dplyr::filter(length >= 0.2 * qlen)

notechis_table <- tibble(qseqid = names(table(blastn_notechis$qseqid)), notechis_n = table(blastn_notechis$qseqid))

print("Notechis scutatus")

joined_counts <- left_join(aipysurus_table, emydocephalus_table) %>%
  left_join(hydrophis_table) %>%
  left_join(notechis_table)

joined_counts %>%
  filter(is.na(notechis_n))

