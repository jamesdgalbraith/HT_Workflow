library(tidyverse)
library(plyranges)
library(BSgenome)

codes <- read_tsv("~/Databases/localrpsb/cddid.tbl3", col_names = c("sseqid", "db_no", "code"))

# Filtering out 
Rex1_seq <- Biostrings::readDNAStringSet(filepath = "TreeBuilding/seqs/All_Rex1.fasta")

Rex1_rps <- read_tsv(file = "TreeBuilding/seqs/All_Rex1_rps.out",
                     col_names = c("seqnames", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs"))  %>%
  dplyr::mutate(strand = case_when(qend > qstart ~ "+", qend < qstart ~ "-"), sseqid = as.double(gsub("gnl\\|CDD\\|", "", sseqid))) %>%
  arrange(seqnames, qstart) %>%
  dplyr::filter((send - sstart + 1) / slen >= 0.5) %>%
  filter(qlen > 800) %>%
  inner_join(codes)

Rex1_RT <- Rex1_rps %>%
  filter(code %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_1", "TERT"))
Rex1_EN <- Rex1_rps %>%
  filter(code %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"))

Rex1_ranges <- Rex1_rps %>%
  dplyr::filter(seqnames %in% Rex1_EN$seqnames, seqnames %in% Rex1_RT$seqnames) %>%
  dplyr::select(seqnames, qlen) %>%
  dplyr::mutate(start = 1) %>%
  dplyr::rename(end = qlen) %>%
  base::unique() %>%
  plyranges::as_granges()

Rex1_filtered_seq <- Biostrings::getSeq(Rex1_seq, Rex1_ranges)

names(Rex1_filtered_seq) <- seqnames(Rex1_ranges)

Biostrings::writeXStringSet(x = Rex1_filtered_seq, filepath = "TreeBuilding/seqs/All_Rex1_filtered.fasta")

RTE_seq <- Biostrings::readDNAStringSet(filepath = "TreeBuilding/seqs/All_RTE.fasta")

RTE_rps <- read_tsv(file = "TreeBuilding/seqs/All_RTE_rps.out",
                col_names = c("seqnames", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs"))  %>%
  dplyr::mutate(strand = case_when(qend > qstart ~ "+", qend < qstart ~ "-"), sseqid = as.double(gsub("gnl\\|CDD\\|", "", sseqid))) %>%
  arrange(seqnames, qstart) %>%
  dplyr::filter((send - sstart + 1) / slen >= 0.5) %>%
  filter(qlen > 800) %>%
  inner_join(codes)

RTE_RT <- RTE_rps %>%
  filter(code %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_1", "TERT"))
RTE_EN <- RTE_rps %>%
  filter(code %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"))

RTE_ranges <- RTE_rps %>%
  filter(seqnames %in% RTE_RT$seqnames, seqnames %in% RTE_EN$seqnames) %>%
  dplyr::select(seqnames, qlen) %>%
  dplyr::mutate(start = 1) %>%
  dplyr::rename(end = qlen) %>%
  base::unique() %>%
  plyranges::as_granges()

RTE_filtered_seq <- Biostrings::getSeq(RTE_seq, RTE_ranges)

names(RTE_filtered_seq) <- seqnames(RTE_ranges)

Biostrings::writeXStringSet(x = RTE_filtered_seq, filepath = "TreeBuilding/seqs/All_RTE_filtered.fasta")
