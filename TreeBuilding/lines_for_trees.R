library(tidyverse)
library(plyranges)
library(BSgenome)

# Filtering out 
Rex1 <- read_tsv(file = "trees/seqs/All_Rex1_rpst.out",
                 col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen")) %>%
  filter((send - sstart + 1) >= slen * 0.8) %>%
  mutate(sseqid = gsub("gnl\\|CDD\\|", "", sseqid))

Rex1_EN <- Rex1 %>%
  filter(sseqid %in% c(197310, 335306, 339261, 197306, 197311))

Rex1_RT <- Rex1 %>%
  filter(sseqid %in% c(238827, 333820, 238828, 275209))

Rex1_ranges <- Rex1 %>%
  dplyr::filter(qseqid %in% Rex1_EN$qseqid, qseqid %in% Rex1_RT$qseqid) %>%
  dplyr::mutate(seqnames = qseqid, start = 1, end = qlen) %>%
  dplyr::select(seqnames, start, end) %>%
  base::unique() %>%
  plyranges::as_granges()

Rex1_seq <- Biostrings::readDNAStringSet(filepath = "trees/seqs/All_Rex1.fasta")

names(Rex1_seq) <- sub(" .*", "", names(Rex1_seq))

Rex1_filtered_seq <- Biostrings::getSeq(Rex1_seq, Rex1_ranges)

names(Rex1_filtered_seq) <- seqnames(Rex1_ranges)

Biostrings::writeXStringSet(x = Rex1_filtered_seq, filepath = "trees/seqs/All_Rex1_filtered.fasta")

RTE <- read_tsv(file = "trees/seqs/All_RTE_rpst.out",
                col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen")) %>%
  filter((send - sstart + 1) >= slen * 0.8) %>%
  mutate(sseqid = gsub("gnl\\|CDD\\|", "", sseqid))

RTE_EN <- RTE %>%
  filter(sseqid %in% c(197310, 335306, 339261, 197306, 197311))

RTE_RT <- RTE %>%
  filter(sseqid %in% c(238827, 333820, 238828, 275209))

RTE_ranges <- RTE %>%
  dplyr::filter(qseqid %in% RTE_EN$qseqid, qseqid %in% RTE_RT$qseqid) %>%
  dplyr::mutate(seqnames = qseqid, start = 1, end = qlen) %>%
  dplyr::select(seqnames, start, end) %>%
  base::unique() %>%
  plyranges::as_granges()

RTE_seq <- Biostrings::readDNAStringSet(filepath = "trees/seqs/All_RTE.fasta")

names(RTE_seq) <- sub(" .*", "", names(RTE_seq))

RTE_filtered_seq <- Biostrings::getSeq(RTE_seq, RTE_ranges)

names(RTE_filtered_seq) <- seqnames(RTE_ranges)

Biostrings::writeXStringSet(x = RTE_filtered_seq, filepath = "trees/seqs/All_RTE_filtered.fasta")