library(tidyverse)
library(GenomicRanges)
library(BSgenome)
library(plyranges)

### NEED TO REPATH PATHS ###

# set species and genome names
species_name <- "Aipysurus_laevis"

# make directory
if(!dir.exists(paste0("CARP/", species_name, "/long_seq"))){dir.create(paste0("CARP/", species_name, "/long_seq"))}

# read in codes
codes <- read_tsv("~/Databases/localrpsb/cddid.tbl3", col_names = c("sseqid", "db_no", "code"))

# read in rps
rps_out <- read_tsv(paste0("CARP/", species_name, "/", species_name, "_rps.out"),
                    col_names = c("seqnames", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs")) %>%
  dplyr::mutate(strand = case_when(qend > qstart ~ "+", qend < qstart ~ "-"), sseqid = as.double(gsub("gnl\\|CDD\\|", "", sseqid))) %>%
  arrange(seqnames, qstart) %>%
  dplyr::filter((send - sstart + 1) / slen >= 0.5) %>%
  filter(qlen > 800) %>%
  inner_join(codes)

# read in CARP
carp_seq <- Biostrings::readDNAStringSet(paste0("CARP/", species_name, "/", species_name, "_Denovo_TE_Library.fasta"))
names(carp_seq) <- sub(" .*", "", names(carp_seq))
carp_tibble <- tibble(seqnames = names(carp_seq), end = width(carp_seq))

long_carp_ranges <- carp_tibble %>%
  dplyr::filter(end >= 800) %>%
  mutate(start = 1) %>%
  plyranges::as_granges()

long_carp_seq <- getSeq(carp_seq, long_carp_ranges)
names(long_carp_seq) <- seqnames(long_carp_ranges)
Biostrings::writeXStringSet(long_carp_seq, paste0("CARP/", species_name, "/long_seq/long_carp_seqs.fa"))
system(paste0("makeblastdb -dbtype nucl -in CARP/", species_name, "/long_seq/long_carp_seqs.fa -out CARP/", species_name, "/long_seq/long_carp_seqs.fa"))

# LINEs
RT <- rps_out %>%
  filter(code %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_1", "TERT"))
EN <- rps_out %>%
  filter(code %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"))

if(nrow(RT) > 0 & nrow(EN) > 0){
  
  potential_LINEs <- rps_out %>%
    filter(seqnames %in% RT$seqnames, seqnames %in% EN$seqnames, qlen < 10000) %>%
    select(seqnames, qlen, strand) %>%
    mutate(start = 1) %>%
    dplyr::rename(end = qlen) %>%
    base::unique() %>%
    as_granges()
  
  potential_LINEs_seq <- getSeq(carp_seq, potential_LINEs)
  names(potential_LINEs_seq) <- seqnames(potential_LINEs)
  writeXStringSet(potential_LINEs_seq, paste0("CARP/", species_name, "/long_seq/potential_LINEs.fa"))
  system(paste0("makeblastdb -dbtype nucl -in CARP/", species_name, "/long_seq/potential_LINEs.fa -out CARP/", species_name, "/long_seq/temp"))
  
  potential_LINEs_blast <- read_tsv(system(paste0("blastn -dust yes -num_threads 12 -db CARP/", species_name, "/long_seq/temp -query CARP/", species_name, "/long_seq/potential_LINEs.fa -outfmt \"6 qseqid sseqid qstart qend sstart send pident qcovs bitscore length mismatch evalue qlen slen\" -task blastn"), intern = TRUE), col_names = c("qseqid", "seqnames", "qstart", "qend", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "qlen", "slen"))
  
  redundant_potential_LINEs <- potential_LINEs_blast %>%
    dplyr::filter(seqnames != qseqid, pident >= 94, qlen < slen, qcovs > 50) %>%
    dplyr::select(seqnames) %>%
    base::unique()
  
  potential_LINEs <- potential_LINEs %>%
    filter(!seqnames %in% redundant_potential_LINEs$seqnames)
  potential_LINEs_seq <- getSeq(carp_seq, potential_LINEs)
  names(potential_LINEs_seq) <- seqnames(potential_LINEs)
  writeXStringSet(potential_LINEs_seq, paste0("CARP/", species_name, "/long_seq/potential_LINEs.fa"))
  
}

# Penelopes
GIY_YIG <- rps_out %>%
  filter(code %in% c("GIY-YIG_PLEs"))

if(nrow(GIY_YIG) > 0){
  
  potential_Penelopes <- rps_out %>%
    filter(seqnames %in% GIY_YIG$seqnames, qlen < 7000) %>%
    select(seqnames, qlen) %>%
    mutate(start = 1) %>%
    dplyr::rename(end = qlen) %>%
    base::unique() %>%
    as_granges()
  
  potential_Penelopes_seq <- getSeq(carp_seq, potential_Penelopes)
  names(potential_Penelopes_seq) <- seqnames(potential_Penelopes)
  writeXStringSet(potential_Penelopes_seq, paste0("CARP/", species_name, "/long_seq/potential_Penelopes.fa"))
  system(paste0("makeblastdb -dbtype nucl -in CARP/", species_name, "/long_seq/potential_Penelopes.fa -out CARP/", species_name, "/long_seq/temp"))
  
  potential_Penelopes_blast <- read_tsv(system(paste0("blastn -dust yes -num_threads 12 -db CARP/", species_name, "/long_seq/temp -query CARP/", species_name, "/long_seq/potential_Penelopes.fa -outfmt \"6 sseqid qseqid qstart qend sstart send pident qcovs bitscore length mismatch evalue qlen slen\" -task blastn"), intern = TRUE), col_names = c("seqnames", "qseqid",  "qstart", "qend", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "qlen", "slen"))
  
  redundant_potential_Penelopes <- potential_Penelopes_blast %>%
    dplyr::filter(seqnames != qseqid, pident >= 94, qlen < slen, qcovs > 75) %>%
    dplyr::select(seqnames) %>%
    base::unique()
  
  potential_Penelopes <- potential_Penelopes %>%
    filter(!seqnames %in% redundant_potential_Penelopes$seqnames)
  potential_Penelopes_seq <- getSeq(carp_seq, potential_Penelopes)
  names(potential_Penelopes_seq) <- seqnames(potential_Penelopes)
  writeXStringSet(potential_Penelopes_seq, paste0("CARP/", species_name, "/long_seq/potential_Penelopes.fa"))
  
}

# PIF/Harbingers
DDE_Tnp_1_4 <- rps_out %>% filter(code %in% c("DDE_Tnp_1", "DDE_Tnp_4", "Plant_tran"))
MADF <- rps_out %>% filter(code %in% c("MADF_DNA_bdg", "MADF", "GT1", "SANT"))

if(nrow(DDE_Tnp_1_4) > 0 & nrow(MADF) > 0){
  
  potential_Harbinger <- rps_out %>%
    filter(seqnames %in% DDE_Tnp_1_4$seqnames, seqnames %in% MADF$seqnames, qlen < 8000) %>%
    select(seqnames, qlen) %>%
    mutate(start = 1) %>%
    dplyr::rename(end = qlen) %>%
    base::unique() %>%
    as_granges()
  
  temp <- rps_out %>%
    filter(seqnames %in% seqnames(potential_Harbinger)) %>%
    select(code)
  
  table(temp$code)
  
  potential_Harbinger_seq <- getSeq(carp_seq, potential_Harbinger)
  names(potential_Harbinger_seq) <- seqnames(potential_Harbinger)
  writeXStringSet(potential_Harbinger_seq, paste0("CARP/", species_name, "/long_seq/potential_Harbinger.fa"))
  system(paste0("makeblastdb -dbtype nucl -in CARP/", species_name, "/long_seq/potential_Harbinger.fa -out CARP/", species_name, "/long_seq/temp"))
  
  potential_Harbinger_blast <- read_tsv(system(paste0("blastn -dust yes -num_threads 12 -db CARP/", species_name, "/long_seq/temp -query CARP/", species_name, "/long_seq/potential_Harbinger.fa -outfmt \"6 qseqid sseqid qstart qend sstart send pident qcovs bitscore length mismatch evalue qlen slen\" -task blastn"), intern = TRUE), col_names = c("qseqid", "seqnames", "qstart", "qend", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "qlen", "slen"))
  
  redundant_potential_Harbinger <- potential_Harbinger_blast %>%
    dplyr::filter(seqnames != qseqid, pident >= 94, qlen < slen, qcovs > 50) %>%
    dplyr::select(seqnames) %>%
    base::unique()
  
  potential_Harbinger <- potential_Harbinger %>%
    filter(!seqnames %in% redundant_potential_Harbinger$seqnames)
  potential_Harbinger_seq <- getSeq(carp_seq, potential_Harbinger)
  names(potential_Harbinger_seq) <- seqnames(potential_Harbinger)
  writeXStringSet(potential_Harbinger_seq, paste0("CARP/", species_name, "/long_seq/potential_Harbinger.fa"))
}

# TcMariners
Tc1_DDE_HTH <- rps_out %>% filter(code %in% c("DDE_1", "DDE_3", "HTH_20", "HTH_23", "HTH_24", "HTH_28", "HTH_32",
                                              "HTH_5", "HTH_7", "HTH_AsnC-type", "HTH_Tnp_ISL3", "HTH_Tnp_Tc3_2"))

if(nrow(Tc1_DDE_HTH) > 0){
  
  potential_TcMar <- rps_out %>%
    filter(seqnames %in% Tc1_DDE_HTH$seqnames, qlen < 5000) %>%
    select(seqnames, qlen, strand) %>%
    mutate(start = 1) %>%
    dplyr::rename(end = qlen) %>%
    base::unique() %>%
    as_granges()
  
  potential_TcMar_seq <- getSeq(carp_seq, potential_TcMar)
  names(potential_TcMar_seq) <- seqnames(potential_TcMar)
  writeXStringSet(potential_TcMar_seq, paste0("CARP/", species_name, "/long_seq/potential_TcMar.fa"))
  system(paste0("makeblastdb -dbtype nucl -in CARP/", species_name, "/long_seq/potential_TcMar.fa -out CARP/", species_name, "/long_seq/temp"))
  
  potential_TcMar_blast <- read_tsv(system(paste0("blastn -dust yes -num_threads 12 -db CARP/", species_name, "/long_seq/temp -query CARP/", species_name, "/long_seq/potential_TcMar.fa -outfmt \"6 sseqid qseqid qstart qend sstart send pident qcovs bitscore length mismatch evalue qlen slen\" -task blastn"), intern = TRUE), col_names = c("seqnames", "qseqid",  "qstart", "qend", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "qlen", "slen"))
  
  redundant_potential_TcMar <- potential_TcMar_blast %>%
    dplyr::filter(seqnames != qseqid, pident >= 94, qlen < slen, qcovs > 75) %>%
    dplyr::select(seqnames) %>%
    base::unique()
  
  potential_TcMar <- potential_TcMar %>%
    filter(!seqnames %in% redundant_potential_TcMar$seqnames)
  potential_TcMar_seq <- getSeq(carp_seq, potential_TcMar)
  names(potential_TcMar_seq) <- seqnames(potential_TcMar)
  writeXStringSet(potential_TcMar_seq, paste0("CARP/", species_name, "/long_seq/potential_TcMar.fa"))
  
}

# hATs
Dimer_Tnp_hAT <- rps_out %>% filter(code == "Dimer_Tnp_hAT")

if(nrow(Dimer_Tnp_hAT) > 0){
  
  potential_hAT <- Dimer_Tnp_hAT %>%
    select(seqnames, qlen, strand) %>%
    mutate(start = 1) %>%
    dplyr::rename(end = qlen) %>%
    base::unique() %>%
    as_granges()
  
  potential_hAT_seq <- getSeq(carp_seq, potential_hAT)
  names(potential_hAT_seq) <- seqnames(potential_hAT)
  writeXStringSet(potential_hAT_seq, paste0("CARP/", species_name, "/long_seq/potential_hAT.fa"))
  system(paste0("makeblastdb -dbtype nucl -in CARP/", species_name, "/long_seq/potential_hAT.fa -out CARP/", species_name, "/long_seq/temp"))
  
  potential_hAT_blast <- read_tsv(system(paste0("blastn -dust yes -num_threads 12 -db CARP/", species_name, "/long_seq/temp -query CARP/", species_name, "/long_seq/potential_hAT.fa -outfmt \"6 sseqid qseqid qstart qend sstart send pident qcovs bitscore length mismatch evalue qlen slen\" -task blastn"), intern = TRUE), col_names = c("seqnames", "qseqid",  "qstart", "qend", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "qlen", "slen"))
  
  redundant_potential_hAT <- potential_hAT_blast %>%
    dplyr::filter(seqnames != qseqid, pident >= 94, qlen < slen, qcovs > 75) %>%
    dplyr::select(seqnames) %>%
    base::unique()
  
  potential_hAT <- potential_hAT %>%
    filter(!seqnames %in% redundant_potential_hAT$seqnames)
  potential_hAT_seq <- getSeq(carp_seq, potential_hAT)
  names(potential_hAT_seq) <- seqnames(potential_hAT)
  writeXStringSet(potential_hAT_seq, paste0("CARP/", species_name, "/long_seq/potential_hAT.fa"))
}

# ERVS
ERV_env <- rps_out %>%
  filter(db_no %in% c("pfam00429"))

if(nrow(ERV_env) > 0){
  
  potential_ERV_elements <- rps_out %>%
    filter(seqnames %in% ERV_env$seqnames) %>%
    dplyr::select(seqnames, qlen, strand) %>%
    dplyr::rename(end = qlen) %>%
    dplyr::mutate(start = 1) %>%
    base::unique() %>%
    plyranges::as_granges()
  
  potential_ERV_elements_seq <- getSeq(carp_seq, potential_ERV_elements)
  names(potential_ERV_elements_seq) <- seqnames(potential_ERV_elements)
  writeXStringSet(potential_ERV_elements_seq, paste0("CARP/", species_name, "/long_seq/potential_ERV_elements.fa"))
  system(paste0("makeblastdb -dbtype nucl -in CARP/", species_name, "/long_seq/potential_ERV_elements.fa -out CARP/", species_name, "/long_seq/temp"))
  
  potential_ERV_elements_blast <- read_tsv(system(paste0("blastn -dust yes -num_threads 12 -db CARP/", species_name, "/long_seq/temp -query CARP/", species_name, "/long_seq/potential_ERV_elements.fa -outfmt \"6 qseqid sseqid qstart qend sstart send pident qcovs bitscore length mismatch evalue qlen slen\" -task blastn"), intern = TRUE), col_names = c("qseqid", "seqnames", "qstart", "qend", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "qlen", "slen"))
  
  redundant_ERV_elements <- potential_ERV_elements_blast %>%
    dplyr::filter(seqnames != qseqid, pident >= 94, qlen < slen, qcovs > 75) %>%
    dplyr::select(qseqid, qlen) %>%
    base::unique()
  
  potential_ERV_elements <- potential_ERV_elements %>%
    filter(!seqnames %in% redundant_ERV_elements$qseqid)
  potential_ERV_elements_seq <- getSeq(carp_seq, potential_ERV_elements)
  names(potential_ERV_elements_seq) <- seqnames(potential_ERV_elements)
  writeXStringSet(potential_ERV_elements_seq, paste0("CARP/", species_name, "/long_seq/potential_ERV_elements.fa"))
  
  potential_extra_ERV_elements <- read_tsv(system(paste0("blastn -dust yes -num_threads 12 -db CARP/", species_name, "/long_seq/long_carp_seqs.fa -query CARP/", species_name, "/long_seq/potential_ERV_elements.fa -outfmt \"6 qseqid sseqid qstart qend sstart send pident qcovs bitscore length mismatch evalue qlen slen\" -task blastn"), intern = TRUE), col_names = c("qseqid", "seqnames", "qstart", "qend", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "qlen", "slen")) %>%
    filter(length > 200) %>%
    select(seqnames, slen) %>%
    base::unique()
  
}

# LTR elements
RT_LTR <- c("RT_DIRS1", "RT_G2_intron", "RT_like", "RT_LTR", "RT_pepA17", "RT_Rtv", "RT_ZFREV_like")

RNase_HI_RT <- rps_out %>%
  filter(grepl("RNase_HI_RT", code))

retropepsin <- rps_out %>%
  filter(qlen > 2000, qlen < 10000) %>%
  filter(code %in% c("retropepsin_like", "retropepsin_like_bacteria", "retropepsin_like_LTR_1", "retropepsin_like_LTR_2", "HIV_retropepsin_like"))

Integrase_H2C2 <- rps_out %>%
  filter(code %in% c("Integrase_H2C2"))

gag  <- rps_out %>%
  filter(code %in% c("gag-asp_proteas", "gag_pre-integrs", "Retrotrans_gag", "Retrotran_gag_2", "Retrotran_gag_3"))

rve  <- rps_out %>%
  filter(code %in% c("rve"))

LTR_temp1 <- rps_out %>%
  filter(qlen > 2000, qlen < 10000) %>%
  filter(seqnames %in% RNase_HI_RT$seqnames, seqnames %in% LTR_RT$seqnames, seqnames %in% gag$seqnames) %>%
  select(seqnames) %>%
  base::unique()

LTR_temp2 <- rps_out %>%
  filter(qlen > 2000, qlen < 10000) %>%
  filter(seqnames %in% rve$seqnames, seqnames %in% Integrase_H2C2$seqnames, seqnames %in% RNase_HI_RT$seqnames) %>%
  select(seqnames) %>%
  base::unique()

LTR_temp3 <- rps_out %>%
  filter(qlen > 2000, qlen < 10000) %>%
  filter(seqnames %in% rve$seqnames & (seqnames %in% Integrase_H2C2$seqnames | seqnames %in% gag$seqnames)) %>%
  select(seqnames) %>%
  base::unique()

if(nrow(LTR_temp1) > 0 | nrow(LTR_temp2) > 0 | nrow(LTR_temp3) > 0 | nrow(retropepsin) > 0){
  
  potential_LTR_elements <- rps_out %>%
    filter(qlen > 2000, qlen < 10000) %>%
    filter(seqnames %in% LTR_temp1$seqnames | seqnames %in% LTR_temp2$seqnames | seqnames %in% LTR_temp3$seqnames | seqnames %in% retropepsin$seqnames) %>%
    base::unique() %>%
    dplyr::select(seqnames, qlen, strand) %>%
    dplyr::rename(end = qlen) %>%
    dplyr::mutate(start = 1) %>%
    base::unique() %>%
    plyranges::as_granges()
  
  if(nrow(potential_extra_ERV_elements > 0)){
    potential_LTR_elements <- potential_LTR_elements %>%
      filter(!seqnames %in% potential_extra_ERV_elements$seqnames)
  }
  
  potential_LTR_elements %>%
    as_tibble() %>%
    arrange(-width)
  
  potential_LTR_elements_seq <- getSeq(carp_seq, potential_LTR_elements)
  names(potential_LTR_elements_seq) <- seqnames(potential_LTR_elements)
  writeXStringSet(potential_LTR_elements_seq, paste0("CARP/", species_name, "/long_seq/potential_LTR_elements.fa"))
  system(paste0("makeblastdb -dbtype nucl -in CARP/", species_name, "/long_seq/potential_LTR_elements.fa -out CARP/", species_name, "/long_seq/temp"))
  
  potential_LTR_elements_blast <- read_tsv(system(paste0("blastn -dust yes -num_threads 12 -db CARP/", species_name, "/long_seq/temp -query CARP/", species_name, "/long_seq/potential_LTR_elements.fa -outfmt \"6 sseqid qseqid qstart qend sstart send pident qcovs bitscore length mismatch evalue qlen slen\" -task blastn"), intern = TRUE), col_names = c("seqnames", "qseqid",  "qstart", "qend", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "qlen", "slen"))
  
  redundant_LTR_elements <- potential_LTR_elements_blast %>%
    dplyr::filter(seqnames != qseqid, pident >= 94, qlen < slen, qcovs > 75) %>%
    dplyr::select(seqnames) %>%
    base::unique()
  
  potential_LTR_elements <- potential_LTR_elements %>%
    filter(!seqnames %in% redundant_LTR_elements$seqnames)
  potential_LTR_elements_seq <- getSeq(carp_seq, potential_LTR_elements)
  names(potential_LTR_elements_seq) <- seqnames(potential_LTR_elements)
  writeXStringSet(potential_LTR_elements_seq, paste0("CARP/", species_name, "/long_seq/potential_LTR_elements.fa"))
}

potentially_classified_seq <- c(potential_ERV_elements_seq, potential_LTR_elements_seq, potential_TcMar_seq,
                                potential_hAT_seq, potential_LINEs_seq, potential_Penelopes_seq)

writeXStringSet(potentially_classified_seq, paste0("CARP/", species_name, "/long_seq/potentially_classified_elements.fa"))



# potentially_classified_blast <- read_tsv(system(paste0("blastn -dust yes -num_threads 12 -db CARP/", species_name, "/long_seq/long_carp_seqs.fa -query CARP/", species_name, "/long_seq/potentially_classified_elements.fa -outfmt \"6 qseqid sseqid qstart qend sstart send pident qcovs bitscore length mismatch evalue qlen slen\" -task blastn"), intern = TRUE), col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "qlen", "slen"))
# 
# potential_LTR_elements %>%
#   as_tibble() %>%
#   arrange(width) %>%
#   as_granges()