# script to identify intact, poorly annotated LINEs in CARP output

library(tidyverse)
library(plyranges)
library(BSgenome)

# read in raw data
genome_fai <- read_tsv(paste0("~/Analysis/Genomes/Aipysurus_laevis/assembly_20171114.fasta.fai"), col_names = c("seqnames", "length", 1:3)) %>%
  select(-c(3:5))

igor <- read_tsv("raw_data/Aipysurus_laevis.igor.gff", col_names = c("seqnames", "X2", "X3", "start", "stop", "X6", "strand", "X8", "Family")) %>%
  select(seqnames, start, stop, Family) %>%
  dplyr::mutate(Family = as.double(gsub("Family ", "", Family)))

consensus <- read_tsv("raw_data/Aipysurus_laevis_Denovo_TE_Library.fasta.fai", col_names = c("seqnames", "end", 1:3)) %>%
  dplyr::mutate(start = 1) %>%
  dplyr::select(seqnames, start, end)
    
# filter for long non repeat and not unknown repeats and short sequences
candidate_names <- consensus %>%
  dplyr::filter(end >= 2300, !grepl("Unclassified", seqnames), !grepl(":", seqnames)) %>%
  dplyr::mutate(seqnames = gsub("#.*", "", seqnames))

# create tibble of above granges object
candidates <- candidate_names  %>%
  plyranges::as_granges()

# read in consensus sequences
consensus_seq <- Biostrings::readDNAStringSet(filepath = "raw_data/Aipysurus_laevis_consensus.fasta")

# rename consensus sequences
names(consensus_seq) <- sub(" \\(.*", "", names(consensus_seq))

# get candidate sequences
candidate_seqs <- Biostrings::getSeq(consensus_seq, candidates)

# name candidate sequences
names(candidate_seqs) <- candidate_names$seqnames

# write candidate sequences to files
Biostrings::writeXStringSet(x = candidate_seqs, filepath = "raw_data/candidate.fa")

# run rpstblastn
system("rpstblastn -db ~/Databases/localrpsb/db/Cdd -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\" -evalue 0.01 -query raw_data/candidate.fa -out raw_data/candidate_cdd.out")

# read in rpstblastn data, rename and filter out hits less that 80% of domain present
rpstblastn_out <- read_tsv("raw_data/candidate_cdd.out", col_names = c("seqnames", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen")) %>%
  dplyr::mutate(strand = case_when(qend > qstart ~ "+", qend < qstart ~ "-"), sseqid = gsub("gnl\\|CDD\\|", "", sseqid)) %>%
  arrange(seqnames, qstart) %>%
  dplyr::filter((send - sstart + 1) / slen >= 0.8)

# determine intact LINEs (sequences containing both RT and EN) and convert to Granges
rpstblastn_en <- rpstblastn_out %>%
  filter(sseqid %in% c(308788, 316995, 197310, 197311)) # endonucleases
rpstblastn_rt <- rpstblastn_out %>%
  dplyr::filter(sseqid %in% c(238827, 238828, 306564)) # reverse transcriptases
LINEs <- rpstblastn_out %>%
  dplyr::filter(seqnames %in% rpstblastn_en$seqnames & seqnames %in% rpstblastn_rt$seqnames) %>%
  select(seqnames, strand) %>%
  base::unique() %>%
  inner_join(candidate_names) %>%
  plyranges::as_granges()

# create names for intact LINEs
LINE_names <- LINEs %>%
  as_tibble() %>%
  mutate(seqnames = as.character(seqnames))

# extract and name intact LINEs seqs from all seq
LINE_seqs <- Biostrings::getSeq(consensus_seq, LINEs)
names(LINE_seqs) <- LINE_names$seqnames