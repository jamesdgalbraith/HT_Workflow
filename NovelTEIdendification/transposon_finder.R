library(tidyverse)
library(GenomicRanges)
library(BSgenome)
library(plyranges)

### NEED TO REPATH PATHS ###

# set species and genome names
species_name <- "Hydrophis_curtis"
genome_path <- "~/Genomes/Reptiles/Aipysurus_laevis/assembly_20171114.fasta"
repbase_path <- "~/Databases/RepBase/RepBase24.07_classed.fasta"
genome_seq <- readDNAStringSet(genome_path)
names(genome_seq) <- sub(" .*", "", names(genome_seq))
gc()



# make directory
if(!dir.exists(paste0("CARP/", species_name, "/long_seq"))){dir.create(paste0("CARP/", species_name, "/long_seq"), recursive = T)}

# read in CARP
carp_seq <- Biostrings::readDNAStringSet(paste0("~/ElapidRepeats/CARP/", species_name, "/consensus.fasta"))
names(carp_seq) <- sub(" .*", "", names(carp_seq))
carp_tibble <- tibble(seqnames = names(carp_seq), end = width(carp_seq))

# create ranges object of long hits
long_carp_ranges <- carp_tibble %>%
  dplyr::filter(end >= 600) %>%
  mutate(start = 1) %>%
  plyranges::as_granges()

# get seq of long repeats, export
long_carp_seq <- getSeq(carp_seq, long_carp_ranges)
names(long_carp_seq) <- seqnames(long_carp_ranges)

if(!dir.exists(paste0("CARP/", species_name, "/long_seq/"))){dir.create(paste0("CARP/", species_name, "/long_seq/"), recursive = T)}
Biostrings::writeXStringSet(long_carp_seq, paste0("CARP/", species_name, "/long_seq/long_carp_seqs.fasta"))

# manually bundled and ran rpstblastn using rpstblaster.sh (cuts up multifasta, searches for predetermined protein domains, compiles results)

# read in rps
rps_out <- read_tsv(paste0("CARP/", species_name, "/long_seq/long_carp_seqs_rps.tsv"),
                    col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs", "stitle")) %>%
  dplyr::mutate(strand = case_when(qend > qstart ~ "+", qend < qstart ~ "-"),
                qseqid = sub("#.*", "", qseqid), qseqid = sub(":.*", "", qseqid)) %>%
  separate(stitle, into = c("code", "name", "description"), sep = ", ") %>%
  dplyr::filter((send - sstart + 1) / slen >= 0.5)

##### Kobolok #####
Kolobok_ranges <- rps_out %>%
  dplyr::filter(name %in% c("THAP", "DM3")) %>%
  dplyr::mutate(seqnames = qseqid, start = 1, end = qlen) %>%
  plyranges::as_granges() %>%
  IRanges::reduce()

if( length(Kolobok_ranges) > 1){
  
  Kolobok_seq <- getSeq(carp_seq, Kolobok_ranges)
  names(Kolobok_seq) <- paste0(seqnames(Kolobok_ranges), "#Kolobok")
  writeXStringSet(Kolobok_seq, paste0("CARP/", species_name, "/long_seq/potential_Kolobok.fa"))
  
  Kolobok_genome_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_Kolobok.fa -db ", genome_path,
                                             " -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12"), intern = T),
                               col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
    filter(length >= 0.5 * qlen)
  
  Kolobok_counts <- tibble(qseqid = names(table(Kolobok_genome_blast$qseqid)), n = as.integer(table(Kolobok_genome_blast$qseqid))) %>%
    filter(n > 2) %>%
    mutate(class = "Kolobok")
  
  Kolobok_carp_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_Kolobok.fa -db CARP/", species_name,
                                           "/long_seq/long_carp_seqs.fasta -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12"),
                                    intern = T), col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
    filter(length >= 0.5 * qlen)
  
  Kolobok_carp_tbl <- Kolobok_carp_blast %>%
    dplyr::select(sseqid, slen) %>%
    base::unique()
  
} else { Kolobok_carp_tbl <- tibble(sseqid = character())}

##### Tc1/Mariner #####
Tc1_HTH <- rps_out %>%
  filter(name %in% c("HTH_Tnp_Tc3_2", "HTH_Tnp_Tc5", "HTH_Tnp_1", "CENPB", "COG3415","HTH_32",
                     "HTH_ARSR", "HTH_28", "HTH_38", "HTH_7", "HTH_48", "HTH_23", "HTH_24"))

Tc1_DDE <- rps_out %>%
  filter(name %in% c("DDE_1", "DDE_3", "Transposase_1", "BrkDBD", "COG3335"))

Tc1_ranges <- rps_out %>%
  filter(qseqid %in% Tc1_DDE$qseqid, qseqid %in% Tc1_HTH$qseqid) %>%
  mutate(seqnames = qseqid, start = 1, end = qlen) %>%
  as_granges() %>%
  reduce()

if(length(Tc1_ranges) > 1){
  
  Tc1_seq <- getSeq(carp_seq, Tc1_ranges)
  names(Tc1_seq) <- seqnames(Tc1_ranges)
  writeXStringSet(Tc1_seq, paste0("CARP/", species_name, "/long_seq/potential_Tc1.fa"))
  
  Tc1_genome_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_Tc1.fa -db ", genome_path,
                                             " -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12"), intern = T),
                                   col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
    filter(length >= 0.5 * qlen)
  
  Tc1_counts <- tibble(qseqid = names(table(Tc1_genome_blast$qseqid)), n = as.integer(table(Tc1_genome_blast$qseqid))) %>%
    filter(n > 2) %>%
    mutate(class = "Tc1_Mariner")
  
  Tc1_carp_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_Tc1.fa -db CARP/", species_name,
                                           "/long_seq/long_carp_seqs.fasta -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12"),
                                    intern = T), col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
    filter(length >= 0.2 * qlen)
  
  Tc1_carp_tbl <- Tc1_carp_blast %>%
    dplyr::select(sseqid, slen) %>%
    base::unique()
  
} else { Tc1_carp_tbl <- tibble(sseqid = character())}

##### hAT #####
hAT_Tnp <- rps_out %>%
  filter(name %in% c("Dimer_Tnp_hAT"))

# hAT_ZnF <- rps_out %>%
#   filter(name %in% c("ZnF_BED", "zf-BED"))

hAT_ranges <- rps_out %>%
  filter(qseqid %in% hAT_Tnp$qseqid) %>%
  # filter(qseqid %in% hAT_ZnF$qseqid) %>%
  mutate(seqnames = qseqid, start = 1, end = qlen) %>%
  as_granges() %>%
  reduce()

if(length(hAT_ranges) > 1){
  
  hAT_seq <- getSeq(carp_seq, hAT_ranges)
  names(hAT_seq) <- seqnames(hAT_ranges)
  writeXStringSet(hAT_seq, paste0("CARP/", species_name, "/long_seq/potential_hAT.fa"))
  
  hAT_genome_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_hAT.fa -db ", genome_path,
                                             " -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12"), intern = T),
                               col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
    filter(length >= 0.5 * qlen)
  
  hAT_counts <- tibble(qseqid = names(table(hAT_genome_blast$qseqid)), n = as.integer(table(hAT_genome_blast$qseqid))) %>%
    filter(n > 2) %>%
    mutate(class = "hAT")
  
  hAT_carp_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_hAT.fa -db CARP/", species_name,
                                           "/long_seq/long_carp_seqs.fasta -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12"),
                                    intern = T), col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
    mutate(coverage = length/qlen) %>%
    filter(length >= 0.2 * qlen)
  
  hAT_carp_tbl <- hAT_carp_blast %>%
    dplyr::select(sseqid, slen) %>%
    base::unique()
  
} else { hAT_carp_tbl = tibble(sseqid = character())}

##### piggyBac #####

piggyBac_ranges <- rps_out %>%
  filter(name %in% c("DDE_Tnp_1_7")) %>%
  mutate(seqnames = qseqid, start = 1, end = qlen) %>%
  as_granges() %>%
  reduce()

if(length(piggyBac_ranges) > 1){
  
  piggyBac_seq <- getSeq(carp_seq, piggyBac_ranges)
  names(piggyBac_seq) <- seqnames(piggyBac_ranges)
  writeXStringSet(piggyBac_seq, paste0("CARP/", species_name, "/long_seq/potential_piggyBac.fa"))
  
  piggyBac_genome_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_piggyBac.fa -db ", genome_path,
                                             " -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12"), intern = T),
                               col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
    filter(length >= 0.5 * qlen)
  
  piggyBac_counts <- tibble(qseqid = names(table(piggyBac_genome_blast$qseqid)), n = as.integer(table(piggyBac_genome_blast$qseqid))) %>%
    filter(n > 2) %>%
    mutate(class = "piggyBac")
  
  piggyBac_carp_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_piggyBac.fa -db CARP/", species_name,
                                           "/long_seq/long_carp_seqs.fasta -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12"),
                                    intern = T), col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
    mutate(coverage = length/qlen) %>%
    filter(length >= 0.2 * qlen)
  
  piggyBac_carp_tbl <- piggyBac_carp_blast %>%
    dplyr::select(sseqid, slen) %>%
    base::unique()
  
} else { piggyBac_carp_tbl = tibble(sseqid = character())}

##### Harbinger #####
Harbinger_GT1 <- rps_out %>%
  filter(name %in% c("Myb_DNA-bind_4", "GT1"))

Harbinger_TnP <- rps_out %>%
  filter(name %in% c("Plant_tran", "DDE_Tnp_1", "DDE_Tnp_4"))

Harbinger_ranges <- rps_out %>%
  filter(qseqid %in% Harbinger_GT1$qseqid) %>%
  filter(qseqid %in% Harbinger_TnP$qseqid) %>%
  mutate(seqnames = qseqid, start = 1, end = qlen) %>%
  as_granges() %>%
  reduce()

if(length(Harbinger_ranges) > 1){
  
  Harbinger_seq <- getSeq(carp_seq, Harbinger_ranges)
  names(Harbinger_seq) <- seqnames(Harbinger_ranges)
  writeXStringSet(Harbinger_seq, paste0("CARP/", species_name, "/long_seq/potential_Harbinger.fa"))
  
  Harbinger_genome_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_Harbinger.fa -db ", genome_path,
                                             " -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12"), intern = T),
                               col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
    filter(length >= 0.5 * qlen)
  
  Harbinger_counts <- tibble(qseqid = names(table(Harbinger_genome_blast$qseqid)), n = as.integer(table(Harbinger_genome_blast$qseqid))) %>%
    filter(n > 2) %>%
    mutate(class = "Harbinger")
  
  Harbinger_carp_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_Harbinger.fa -db CARP/", species_name,
                                           "/long_seq/long_carp_seqs.fasta -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12"),
                                    intern = T), col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
    mutate(coverage = length/qlen) %>%
    filter(length >= 0.2 * qlen)
  
  Harbinger_carp_tbl <- Harbinger_carp_blast %>%
    dplyr::select(sseqid, slen) %>%
    base::unique()
  
} else { Harbinger_carp_tbl = tibble(sseqid = character())}

##### ERV #####
ERV_env <- rps_out %>%
  filter(name %in% c("Ebola_HIV-1-like_HR1-HR2", "Ebola_RSV-like_HR1-HR2", "Ebola-like_HR1-HR2", "HIV-1-like_HR1-HR2", "GP41",
                     "ENVV1-like_HR1-HR2", "HERV-Rb-like_HR1-HR2", "HTLV-1-like_HR1-HR2", "RSV-like_HR1-HR2", "TLV_coat"))

# ERV_RNase <- rps_out %>%
#   filter(name %in% c("PRK08719", "RNase_H", "RNase_H_Dikarya_like", "RNase_H_like", "RNase_HI_bacteria_like", "RNase_HI_eukaryote_like",
#                      "RNase_HI_like RNase_HI_prokaryote_like", "RNase_HI_RT_Bel", "RNase_HI_RT_DIRS1", "Rnase_HI_RT_non_LTR", "rnhA", "RnhA"))
# 
# ERV_RT <- rps_out %>%
#   filter(name %in% c("RT_ZFREV_like", "RT_Rtv", "RT_LTR", "RT_DIRS1", "RVT_1", "RT_like"))
# 
# ERV_Integrase <- rps_out %>%
#   filter(name %in% c("Integrase_H2C2", "zf-H2C2", "zf-H3C2"))
# 
rps_out %>%
  filter(qseqid %in% ERV_env$qseqid) %>%
  select(qseqid, name, qlen)

ERV_ranges <- rps_out %>%
  filter(qseqid %in% ERV_env$qseqid) %>%
  # filter(qseqid %in% ERV_RT$qseqid, qseqid %in% ERV_RNase$qseqid, ERV_Integrase$qseqid) %>%
  mutate(seqnames = qseqid, start = 1, end = qlen) %>%
  as_granges() %>%
  reduce()

if(length(ERV_ranges) > 1){
  
  ERV_seq <- getSeq(carp_seq, ERV_ranges)
  names(ERV_seq) <- seqnames(ERV_ranges)
  writeXStringSet(ERV_seq, paste0("CARP/", species_name, "/long_seq/potential_ERV.fa"))
  
  ERV_genome_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_ERV.fa -db ", genome_path,
                                                   " -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12"), intern = T),
                                     col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
    filter(length >= 0.5 * qlen)
  
  ERV_counts <- tibble(qseqid = names(table(ERV_genome_blast$qseqid)), n = as.integer(table(ERV_genome_blast$qseqid))) %>%
    filter(n > 2) %>%
    mutate(class = "ERV")
  
  ERV_carp_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_ERV.fa -db CARP/", species_name,
                                                 "/long_seq/long_carp_seqs.fasta -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12"),
                                          intern = T), col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
    mutate(coverage = length/qlen) %>%
    filter(length >= 0.2 * qlen)
  
  ERV_carp_tbl <- ERV_carp_blast %>%
    dplyr::select(sseqid, slen) %>%
    base::unique() %>%
    dplyr::arrange(-slen)
  
}

##### Copia #####
Copia_rve <- rps_out %>%
  filter(name %in% c("rve", "rve3"))

Copia_RT <- rps_out %>%
  filter(name %in% c("RVT_2"))

Copia_RNase <- rps_out %>%
  filter(name %in% c("RNase_HI_RT_Ty1"))

Copia_ranges <- rps_out %>%
  filter(qseqid %in% Copia_RT$qseqid, qseqid %in% Copia_RNase$qseqid, qseqid %in% Copia_rve$qseqid) %>%
  mutate(seqnames = qseqid, start = 1, end = qlen) %>%
  as_granges() %>%
  reduce()

if(length(Copia_ranges) > 1){
  
  Copia_seq <- getSeq(carp_seq, Copia_ranges)
  names(Copia_seq) <- seqnames(Copia_ranges)
  writeXStringSet(Copia_seq, paste0("CARP/", species_name, "/long_seq/potential_Copia.fa"))
  
  Copia_genome_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_Copia.fa -db ", genome_path,
                                             " -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12"), intern = T),
                               col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
    filter(length >= 0.5 * qlen)
  
  Copia_counts <- tibble(qseqid = names(table(Copia_genome_blast$qseqid)), n = as.integer(table(Copia_genome_blast$qseqid))) %>%
    filter(n > 2) %>%
    mutate(class = "Copia")
  
  Copia_carp_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_Copia.fa -db CARP/", species_name,
                                           "/long_seq/long_carp_seqs.fasta -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12"),
                                    intern = T), col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
    mutate(coverage = length/qlen) %>%
    filter(length >= 0.2 * qlen)
  
  Copia_carp_tbl <- Copia_carp_blast %>%
    dplyr::select(sseqid, slen) %>%
    base::unique() %>%
    dplyr::arrange(-slen)
  
}

rps_out %>% filter(qseqid == "family081324")

##### LTR #####
LTR_RNase <- rps_out %>%
  filter(name %in% c("RT_RNaseH_2", "RT_RNaseH", "RNase_HI_RT_Ty3", "RNase_HI_RT_DIRS1", "RNase_H_like"))

LTR_RT <- rps_out %>%
  filter(name %in% c("RT_ZFREV_like", "RT_Rtv", "RT_LTR", "RT_DIRS1", "RVT_1", "RT_like"))

LTR_Integrase <- rps_out %>%
  filter(name %in% c("Integrase_H2C2", "zf-H2C2"))

LTR_ranges <- rps_out %>%
  filter(qseqid %in% LTR_RNase$qseqid, qseqid %in% LTR_RT$qseqid, qseqid %in% LTR_Integrase$qseqid) %>%
  filter(!qseqid %in% Copia_carp_tbl$sseqid, !qseqid %in% ERV_carp_tbl$sseqid) %>%
  mutate(seqnames = qseqid, start = 1, end = qlen) %>%
  as_granges() %>%
  reduce()

if(length(LTR_ranges) > 1){
  
  LTR_seq <- getSeq(carp_seq, LTR_ranges)
  names(LTR_seq) <- seqnames(LTR_ranges)
  writeXStringSet(LTR_seq, paste0("CARP/", species_name, "/long_seq/potential_LTR.fa"))
  
  LTR_genome_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_LTR.fa -db ", genome_path,
                                             " -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12"), intern = T),
                               col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
    filter(length >= 0.5 * qlen)
  
  LTR_counts <- tibble(qseqid = names(table(LTR_genome_blast$qseqid)), n = as.integer(table(LTR_genome_blast$qseqid))) %>%
    filter(n > 2) %>%
    mutate(class = "LTR")
  
  LTR_carp_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_LTR.fa -db CARP/", species_name,
                                           "/long_seq/long_carp_seqs.fasta -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12"),
                                    intern = T), col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
    mutate(coverage = length/qlen) %>%
    filter(length >= 0.2 * qlen)
  
  LTR_carp_tbl <- LTR_carp_blast %>%
    dplyr::select(sseqid, slen) %>%
    base::unique() %>%
    dplyr::arrange(slen)
  
}

LTR_carp_blast %>%
  filter(qseqid == "family159136_consensus", sseqid %in% LTR_counts$qseqid)

##### LINEs #####
LINE_RT <- rps_out %>%
  filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_1", "TERT"))
LINE_EN <- rps_out %>%
  filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"))

# select LINEs based on the presence of both RT and EN, and be shorter than 10000bp
RT_EN <- rps_out %>%
  filter(qseqid %in% LINE_RT$qseqid, qseqid %in% LINE_EN$qseqid, qlen < 10000) %>%
  dplyr::rename(seqnames = qseqid) %>%
  select(seqnames, qlen, strand) %>%
  base::unique()

if(nrow(RT_EN) > 0){
  
  # make rabges object of potential LINEs
  potential_LINEs <- RT_EN %>%
    mutate(start = 1) %>%
    dplyr::rename(end = qlen) %>%
    arrange(seqnames) %>%
    as_granges()
  
  # get seq of potential LINEs, name and write to file
  potential_LINEs_seq <- getSeq(carp_seq, potential_LINEs)
  names(potential_LINEs_seq) <- paste0(seqnames(potential_LINEs), "(", strand(potential_LINEs), ")")
  writeXStringSet(potential_LINEs_seq, paste0("CARP/", species_name, "/long_seq/potential_LINEs.fa"))
  
  LINE_genome_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_LINEs.fa -db ", genome_path,
                                             " -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12"), intern = T),
                               col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
    filter(length >= 0.5 * qlen)
  
  LINE_counts <- tibble(qseqid = names(table(LINE_genome_blast$qseqid)), n = as.integer(table(LINE_genome_blast$qseqid))) %>%
    filter(n > 2) %>%
    mutate(class = "LINE")
  
  LINE_carp_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_LINEs.fa -db CARP/", species_name,
                                            "/long_seq/long_carp_seqs.fasta -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12"),
                                     intern = T), col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
    mutate(coverage = length/qlen) %>%
    filter(length >= 0.2 * qlen)
  
  LINE_carp_tbl <- LINE_carp_blast %>%
    dplyr::select(sseqid, slen) %>%
    base::unique() %>%
    dplyr::arrange(slen)
  
  # potential_LINEs_selfblast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_LINEs.fa -subject CARP/", species_name, "/long_seq/potential_LINEs.fa -outfmt \"6 qseqid sseqid length pident qlen slen qcovs\""), intern = T),
  #                                       col_names = c("qseqid", "sseqid", "length", "pident", "qlen", "slen", "qcovs"))
  # 
  # potential_LINEs_selfblast <- potential_LINEs_selfblast %>%
  #   filter(qseqid != sseqid)
  # 
  # redundant_potential_LINEs <- potential_LINEs_selfblast %>%
  #   filter(pident >= 94, qlen<slen, qcovs >=50)
  # 
  # redundant_potential_LINEs <- redundant_potential_LINEs %>%
  #   dplyr::select(qseqid) %>%
  #   base::unique()
  # 
  # potential_LINEs <- RT_EN %>%
  #   filter(!seqnames %in% redundant_potential_LINEs) %>%
  #   mutate(start = 1) %>%
  #   dplyr::rename(end = qlen) %>%
  #   arrange(seqnames) %>%
  #   as_granges()
  # 
  # potential_LINEs_seq <- getSeq(carp_seq, potential_LINEs)
  # names(potential_LINEs_seq) <- paste0(seqnames(potential_LINEs), "(", strand(potential_LINEs), ")")
  # writeXStringSet(potential_LINEs_seq, paste0("CARP/", species_name, "/long_seq/potential_LINEs.fa"))
 
}

if(length(potential_LINEs) > 0){
  
  # search againist repbase for the most similar classified repeat
  repbase_blast_out <- read_tsv(system(paste0("blastn -evalue 0.00002 -word_size 7 -dust yes -num_threads 12 -reward 3 -penalty -4 -gapopen 30 -gapextend 6 -xdrop_ungap 80 -xdrop_gap 130 -xdrop_gap_final 150 -query CARP/", species_name, "/long_seq/potential_LINEs.fa -db ", repbase_path, " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs")) %>%
    group_by(qseqid) %>%
    arrange(-bitscore) %>%
    dplyr::slice(1) %>%
    ungroup()
  
  # search genome for blast
  genome_blast_LINEs <- read_tsv(system(paste0("blastn -dust yes -num_threads 12 -query CARP/", species_name, "/long_seq/potential_LINEs.fa -db ", genome_path, " -outfmt \"6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen qcovs\" -perc_identity 94"), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs"))
  
  new_LINE_names <- repbase_blast_out %>%
    dplyr::select(qseqid, seqnames) %>%
    dplyr::mutate(new_name = paste0(sub("#.*", "#", qseqid), sub("#.*", "", seqnames))) %>%
    dplyr::select(-seqnames)
  
  # select best hits to each blast search
  genome_blast_LINEs_filtered <- genome_blast_LINEs %>%
    group_by(qseqid) %>%
    arrange(-bitscore) %>%
    dplyr::slice(1:20) %>%
    ungroup()

  if(!dir.exists(paste0("CARP/", species_name, "/long_seq/curation/LINEs"))){dir.create(paste0("CARP/", species_name, "/long_seq/curation/LINEs"), recursive = T)}
  
  # individually align each LINE subfamily
  for(i in 1:length(potential_LINEs)){
    
    # select LINE i
    genome_blast_LINEs_filtered_ranges <- genome_blast_LINEs_filtered %>%
      filter(qseqid == as.character(names(potential_LINEs_seq)[i]), length >= qlen/3, pident >= 94) %>%
      dplyr::select(seqnames, sstart, send, slen) %>%
      mutate(strand = ifelse(sstart < send, "+", "-"),
             start = ifelse(sstart < send, sstart - 2000, send - 2000),
             end = ifelse(sstart > send, sstart + 2000, send + 2000),
             start = ifelse(start < 1, 1, start),
             end = ifelse(end > slen, slen, end)) %>%
      as_granges()
    
    if(length(genome_blast_LINEs_filtered_ranges) > 1){
    
    # get seq of LINE i hits
    genome_blast_LINEs_filtered_seq <- getSeq(genome_seq, genome_blast_LINEs_filtered_ranges)
    names(genome_blast_LINEs_filtered_seq) <- paste0(seqnames(genome_blast_LINEs_filtered_ranges), ":", ranges(genome_blast_LINEs_filtered_ranges),
                                                     "(", strand(genome_blast_LINEs_filtered_ranges), ")")
    
    # add query LINE i seq
    genome_blast_LINEs_filtered_seq <- c(potential_LINEs_seq[i], genome_blast_LINEs_filtered_seq)
    
    # write LINE i seqs to file
    writeXStringSet(genome_blast_LINEs_filtered_seq, paste0("CARP/", species_name, "/long_seq/curation/temp.fa"))
    
    # create name for LINE i
    query_repeat_subbed <- paste0(species_name, "_", sub("\\(.\\)", "_", new_LINE_names$new_name[i]))
    
    # align LINE i
    system(paste0("mafft --adjustdirection --localpair --thread 12 CARP/", species_name, "/long_seq/curation/temp.fa > CARP/", species_name, "/long_seq/curation/LINEs/", species_name, "_", query_repeat_subbed, "_aligned.fa"))
    
    }
  }
  
}


##### Penelopes #####
# select GIY_YIG based on RT and GIY_YIG domains
GIY_YIG <- rps_out %>%
  filter(name %in% c("GIY-YIG_PLEs"), qseqid %in% LINE_RT$qseqid, qlen > 2000, qlen < 10000) %>%
  dplyr::select(qseqid, qlen, strand) %>%
  arrange(qseqid) %>%
  dplyr::rename(seqnames = qseqid) %>%
  base::unique()

if(nrow(GIY_YIG) > 0){

  # make rabges object of potential Penelopes
  potential_Penelopes <- GIY_YIG %>%
    mutate(start = 1) %>%
    dplyr::rename(end = qlen) %>%
    arrange(seqnames) %>%
    as_granges()

  # get seq of potential Penelopes, name and write to file
  potential_Penelopes_seq <- getSeq(carp_seq, potential_Penelopes)
  names(potential_Penelopes_seq) <- paste0(seqnames(potential_Penelopes), "(", strand(potential_Penelopes), ")")
  writeXStringSet(potential_Penelopes_seq, paste0("CARP/", species_name, "/long_seq/potential_Penelopes.fa"))

  Penelopes_genome_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_Penelopes.fa -db ", genome_path,
                                              " -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12"), intern = T),
                                col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
    filter(length >= 0.5 * qlen)
  
  Penelopes_counts <- tibble(qseqid = names(table(Penelopes_genome_blast$qseqid)), n = as.integer(table(Penelopes_genome_blast$qseqid))) %>%
    filter(n > 2) %>%
    mutate(class = "Penelope")
  
  Penelopes_carp_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_Penelopes.fa -db CARP/", species_name,
                                            "/long_seq/long_carp_seqs.fasta -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12"),
                                     intern = T), col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
    mutate(coverage = length/qlen) %>%
    filter(length >= 0.2 * qlen)
  
  Penelopes_carp_tbl <- Penelopes_carp_blast %>%
    dplyr::select(sseqid, slen) %>%
    base::unique() %>%
    dplyr::arrange(slen)
}  

compiled_counts <- rbind(hAT_counts, Tc1_counts, ERV_counts, Copia_counts, LTR_counts, LINE_counts, Penelopes_counts) %>%
  mutate(qseqid = sub("\\(.*", "", qseqid))

compiled_counts_ranges <- rps_out %>%
  dplyr::group_by(qseqid) %>%
  dplyr::arrange(-bitscore) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(qseqid, qlen, strand) %>%
  dplyr::inner_join(compiled_counts) %>%
  dplyr::mutate(start = 1) %>%
  dplyr::rename(seqnames = qseqid, end = qlen) %>%
  as_granges()

compiled_counts_seq <- getSeq(carp_seq, compiled_counts_ranges)
names(compiled_counts_seq) <- seqnames(compiled_counts_ranges)
writeXStringSet(compiled_counts_seq, paste0("CARP/", species_name, "/long_seq/potential_compiled.fa"))

compiled_Hydrophis_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_compiled.fa -db ~/Genomes/Reptiles/Hydrophis_curtis/Hcur1.v1.1.fasta -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12 -task dc-megablast"),
                                           intern = T), col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
  mutate(coverage = length/qlen) 

compiled_Hydrophis_tbl <- tibble(qseqid = names(table(compiled_Hydrophis_blast$qseqid)), n = as.integer(table(compiled_Hydrophis_blast$qseqid)))

absent_Hydrophis <- compiled_counts_ranges %>%
  filter(!seqnames(compiled_counts_ranges) %in% compiled_Hydrophis_tbl$qseqid)

absent_Hydrophis_seq <- getSeq(carp_seq, absent_Hydrophis)
names(absent_Hydrophis_seq) <- paste0(seqnames(absent_Hydrophis), "_", absent_Hydrophis$class)
writeXStringSet(absent_Hydrophis_seq, paste0("CARP/", species_name, "/long_seq/absent_Hydrophis.fa"))

compiled_Notechis_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_compiled.fa -db ~/Genomes/Reptiles/Notechis_scutatus/TS10Xv2-PRI.fasta -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12 -task dc-megablast"),
                intern = T), col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
  mutate(coverage = length/qlen) 

compiled_Notechis_tbl <- tibble(qseqid = names(table(compiled_Notechis_blast$qseqid)), n = as.integer(table(compiled_Notechis_blast$qseqid)))

absent_Notechis <- compiled_counts_ranges %>%
  filter(!seqnames(compiled_counts_ranges) %in% compiled_Notechis_tbl$qseqid)

absent_Notechis_seq <- getSeq(carp_seq, absent_Notechis)
names(absent_Notechis_seq) <- paste0(seqnames(absent_Notechis), "_", absent_Notechis$class)
writeXStringSet(absent_Notechis_seq, paste0("CARP/", species_name, "/long_seq/absent_Notechis.fa"))


compiled_Pseudonaja_blast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_compiled.fa -db ~/Genomes/Reptiles/Pseudonaja_textilis/EBS10Xv2-PRI.fasta -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12 -task dc-megablast"),
                                           intern = T), col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
  mutate(coverage = length/qlen)

compiled_Pseudonaja_tbl <- tibble(qseqid = names(table(compiled_Pseudonaja_blast$qseqid)), n = as.integer(table(compiled_Pseudonaja_blast$qseqid)))

absent_Pseudonaja <- compiled_counts_ranges %>%
  filter(!seqnames(compiled_counts_ranges) %in% compiled_Pseudonaja_tbl$qseqid)



# identify and remove potential duplicates
  # potential_Penelopes_selfblast <- read_tsv(system(paste0("blastn -query CARP/", species_name, "/long_seq/potential_Penelopes.fa -subject CARP/", species_name, "/long_seq/potential_Penelopes.fa -outfmt \"6 qseqid sseqid length pident qlen slen qcovs\""), intern = T),
  #                                           col_names = c("qseqid", "sseqid", "length", "pident", "qlen", "slen", "qcovs"))
  # potential_Penelopes_selfblast <- potential_Penelopes_selfblast %>%
  #   filter(qseqid != sseqid)
  # redundant_potential_Penelopes <- potential_Penelopes_selfblast %>%
  #   filter(pident >= 94, qlen<slen, qcovs >=50)
  # redundant_potential_Penelopes <- redundant_potential_Penelopes %>%
  #   dplyr::select(qseqid) %>%
  #   base::unique()
  # potential_Penelopes <- GIY_YIG %>%
  #   filter(!seqnames %in% redundant_potential_Penelopes) %>%
  #   mutate(start = 1) %>%
  #   dplyr::rename(end = qlen) %>%
  #   arrange(seqnames) %>%
  #   as_granges()
  # potential_Penelopes_seq <- getSeq(carp_seq, potential_Penelopes)
  # names(potential_Penelopes_seq) <- paste0(seqnames(potential_Penelopes), "(", strand(potential_Penelopes), ")")
  # writeXStringSet(potential_Penelopes_seq, paste0("CARP/", species_name, "/long_seq/potential_Penelopes.fa"))
  # 
  # # search againist repbase for the most similar classified repeat
  # repbase_blast_out <- read_tsv(system(paste0("blastn -evalue 0.00002 -word_size 7 -dust yes -num_threads 12 -reward 3 -penalty -4 -gapopen 30 -gapextend 6 -xdrop_ungap 80 -xdrop_gap 130 -xdrop_gap_final 150 -query CARP/", species_name, "/long_seq/potential_Penelopes.fa -db ", repbase_path, " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs")) %>%
  #   group_by(qseqid) %>%
  #   arrange(-bitscore) %>%
  #   dplyr::slice(1) %>%
  #   ungroup()
  # 
  # # search genome for blast
  # genome_blast_Penelopes <- read_tsv(system(paste0("blastn -dust yes -num_threads 12 -query CARP/", species_name, "/long_seq/potential_Penelopes.fa -db ", genome_path, " -outfmt \"6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen qcovs\" -perc_identity 94"), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs"))
  # 
  # new_Penelope_names <- repbase_blast_out %>%
  #   dplyr::select(qseqid, seqnames) %>%
  #   dplyr::mutate(new_name = paste0(sub("#.*", "#", qseqid), sub("#.*", "", seqnames))) %>%
  #   dplyr::select(-seqnames)
  # 
  # # select best hits to each blast search
  # genome_blast_Penelopes_filtered <- genome_blast_Penelopes %>%
  #   group_by(qseqid) %>%
  #   arrange(-bitscore) %>%
  #   dplyr::slice(1:20) %>%
  #   ungroup()

  if(!dir.exists(paste0("CARP/", species_name, "/long_seq/curation/Penelopes"))){dir.create(paste0("CARP/", species_name, "/long_seq/curation/Penelopes"), recursive = T)}

  # individually align each Penelope subfamily
  # for(i in 1:length(potential_Penelopes)){

    # select Penelope i
    genome_blast_Penelopes_filtered_ranges <- genome_blast_Penelopes_filtered %>%
      filter(qseqid == as.character(names(potential_Penelopes_seq)[i]), length >= qlen/3) %>%
      dplyr::select(seqnames, sstart, send, slen) %>%
      mutate(strand = ifelse(sstart < send, "+", "-"),
             start = ifelse(sstart < send, sstart - 2000, send - 2000),
             end = ifelse(sstart > send, sstart + 2000, send + 2000),
             start = ifelse(start < 1, 1, start),
             end = ifelse(end > slen, slen, end)) %>%
      as_granges()

    # get seq of Penelope i hits
    genome_blast_Penelopes_filtered_seq <- getSeq(genome_seq, genome_blast_Penelopes_filtered_ranges)
    names(genome_blast_Penelopes_filtered_seq) <- paste0(seqnames(genome_blast_Penelopes_filtered_ranges), ":", ranges(genome_blast_Penelopes_filtered_ranges),
                                                         "(", strand(genome_blast_Penelopes_filtered_ranges), ")")

    # add query Penelope i seq
    genome_blast_Penelopes_filtered_seq <- c(potential_Penelopes_seq[i], genome_blast_Penelopes_filtered_seq)

    # write Penelope i seqs to file
    writeXStringSet(genome_blast_Penelopes_filtered_seq, paste0("CARP/", species_name, "/long_seq/curation/temp.fa"))

    # create name for Penelope i
    query_repeat_subbed <- paste0(species_name, "_", sub("\\(.\\)", "_", new_Penelope_names$new_name[i]))

    # align Penelope i
    system(paste0("mafft --adjustdirection --localpair --thread 12 CARP/", species_name, "/long_seq/curation/temp.fa > CARP/", species_name, "/long_seq/curation/Penelopes/", species_name, "_", query_repeat_subbed, "_aligned.fa"))

  }
  # }


##### Unclassified #####
# yet_unclassified <- rps_out %>%
#   filter(!qseqid %in% c(LTR_carp_tbl$sseqid, Harbinger_carp_tbl$sseqid, ERV_carp_tbl$sseqid, Copia_carp_tbl$sseqid, Tc1_carp_tbl$sseqid,
#                         )) %>%
#   arrange(qseqid, qstart) %>%
#   select(qseqid, name, qstart, qend, qlen)
# 
# yet_unclassified_ranges <- yet_unclassified %>%
#   mutate(seqnames = qseqid, start = 1, end = qlen) %>%
#   as_granges()
# 
# 
# yet_unclassified_Integrase_ranges <- yet_unclassified_ranges %>%
#   filter(name == "Integrase_H2C2")
# 
# yet_unclassified_Integrase_seq <- getSeq(carp_seq, yet_unclassified_Integrase_ranges)
# names(yet_unclassified_Integrase_seq) <- seqnames(yet_unclassified_Integrase_ranges)
# writeXStringSet(yet_unclassified_Integrase_seq, paste0("CARP/", species_name, "/long_seq/yet_unclassified_H2C2.fa"))
# 
# yet_unclassified_repbase_blast <- read_tsv(system(paste0("tblastx -query CARP/", species_name, "/long_seq/yet_unclassified_H2C2.fa -db ~/Databases/RepBase/RepBase24.07_classed.fasta -outfmt \"6 qseqid sseqid pident length qlen slen qcovs sstart send\" -num_threads 12"),
#                                                   intern = T), col_names = c("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs", "sstart", "send")) %>%
#   mutate(coverage = length/qlen)
# 
# yet_unclassified_repbase_blast %>%
#   group_by(sseqid) %>%
#   dplyr::slice(1) %>%
#   ungroup()