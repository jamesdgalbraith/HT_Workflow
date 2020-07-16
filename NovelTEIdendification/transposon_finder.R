library(tidyverse)
library(GenomicRanges)
library(BSgenome)
library(plyranges)

### NEED TO REPATH PATHS ###

# set species and genome names
species_name <- "Aipysurus_laevis"
genome_path <- "~/Genomes/Reptiles/Aipysurus_laevis/kmer_49.pilon_x2.sorted.fasta"
repbase_path <- "~/Databases/RepBase/RepBase24.07_classed.fasta"
genome_seq <- readDNAStringSet(genome_path)
names(genome_seq) <- sub(" .*", "", names(genome_seq))
gc()



# make directory
if(!dir.exists(paste0("CARP/", species_name, "/long_seq"))){dir.create(paste0("CARP/", species_name, "/long_seq"))}

# read in codes
codes <- read_tsv("~/Databases/localrpsb/cddid.tbl3", col_names = c("sseqid", "db_no", "code"))

# read in CARP
carp_seq <- Biostrings::readDNAStringSet(paste0("~/ElapidRepeats/CARP/", species_name, "/", species_name, "_Denovo_TE_Library.fasta"))
names(carp_seq) <- sub(" .*", "", names(carp_seq))
carp_tibble <- tibble(seqnames = names(carp_seq), end = width(carp_seq))

# create ranges object of long hits
long_carp_ranges <- carp_tibble %>%
  dplyr::filter(end >= 800) %>%
  mutate(start = 1) %>%
  plyranges::as_granges()

# get seq of long repeats, export, makedb and selfblast
long_carp_seq <- getSeq(carp_seq, long_carp_ranges)
names(long_carp_seq) <- seqnames(long_carp_ranges)
Biostrings::writeXStringSet(long_carp_seq, paste0("CARP/", species_name, "/long_seq/long_carp_seqs.fa"))
system(paste0("makeblastdb -dbtype nucl -in CARP/", species_name, "/long_seq/long_carp_seqs.fa -out CARP/", species_name, "/long_seq/long_carp_seqs.fa"))
system(paste0("blastn -query CARP/", species_name, "/long_seq/long_carp_seqs.fa -db CARP/", species_name, "/long_seq/long_carp_seqs.fa -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\" -out CARP/", species_name, "/long_seq/long_carp_seqs.out"))

long_carp_self_out <- read_tsv(paste0("CARP/", species_name, "/long_seq/long_carp_seqs.out"), col_names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qlen","slen","qcovs"))

long_carp_self_out <- long_carp_self_out %>%
  filter(qseqid != sseqid)

redundant <- long_carp_self_out %>%
  filter(qlen < slen, pident >= 94, qcovs >= 50) %>%
  dplyr::select(qseqid) %>%
  base::unique()

non_redundant <- long_carp_self_out %>%
  filter(!qseqid %in% redundant$qseqid, !sseqid %in% redundant$qseqid) %>%
  dplyr::select(qseqid) %>%
  base::unique()

non_redundant_long_carp_ranges <- carp_tibble %>%
  filter(seqnames %in% non_redundant$qseqid) %>%
  dplyr::filter(end >= 800) %>%
  mutate(start = 1) %>%
  plyranges::as_granges()

long_carp_seq <- getSeq(carp_seq, non_redundant_long_carp_ranges)
names(long_carp_seq) <- seqnames(non_redundant_long_carp_ranges)
Biostrings::writeXStringSet(long_carp_seq, paste0("CARP/", species_name, "/long_seq/non_redundant_long_carp_seqs.fa"))


# read in rps
system(command = paste0("rpstblastn -query CARP/", species_name, "/long_seq/non_redundant_long_carp_seqs.fa -db ~/Databases/localrpsb/transposons/transposons -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs stitle\" -evalue 0.01 -out CARP/", species_name, "/long_carp_seqs_rps.out -num_threads 12"))

rps_out <- read_tsv(paste0("CARP/", species_name, "/long_carp_seqs_rps.out"),
                    col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs", "stitle")) %>%
  dplyr::mutate(strand = case_when(qend > qstart ~ "+", qend < qstart ~ "-")) %>%
  separate(stitle, into = c("code", "name", "description"), sep = ", ") %>%
  dplyr::filter((send - sstart + 1) / slen >= 0.5)

# LINEs
RT <- rps_out %>%
  filter(name %in% c("RT_like", "RT_nLTR_like", "RVT_1", "RT_G2_intron", "RVT_1", "TERT"))
EN <- rps_out %>%
  filter(name %in% c("EEP", "EEP-2", "Exo_endo_phos", "Exo_endo_phos_2", "L1-EN", "R1-I-EN"))

RT_EN <- rps_out %>%
  filter(qseqid %in% RT$qseqid, qseqid %in% EN$qseqid)

if(nrow(RT_EN) > 0){
  
  potential_LINEs <- RT_EN %>%
    dplyr::rename(seqnames = qseqid) %>%
    select(seqnames, qlen, strand) %>%
    base::unique() %>%
    mutate(start = 1) %>%
    dplyr::rename(end = qlen) %>%
    arrange(seqnames) %>%
    as_granges()
  
  potential_LINEs_seq <- getSeq(carp_seq, potential_LINEs)
  names(potential_LINEs_seq) <- seqnames(potential_LINEs)
  writeXStringSet(potential_LINEs_seq, paste0("CARP/", species_name, "/long_seq/potential_LINEs.fa"))
  
}

# Penelopes
GIY_YIG <- rps_out %>%
  filter(name %in% c("GIY-YIG_PLEs"), qseqid %in% RT$qseqid)

if(nrow(GIY_YIG) > 0){
  
  potential_Penelopes <- GIY_YIG %>%
    dplyr::rename(seqnames = qseqid) %>%
    select(seqnames, qlen) %>%
    mutate(start = 1) %>%
    dplyr::rename(end = qlen) %>%
    base::unique() %>%
    as_granges()
  
  potential_Penelopes_seq <- getSeq(carp_seq, potential_Penelopes)
  names(potential_Penelopes_seq) <- seqnames(potential_Penelopes)
  writeXStringSet(potential_Penelopes_seq, paste0("CARP/", species_name, "/long_seq/potential_Penelopes.fa"))
  
}

# PIF/Harbingers
DDE_Tnp_1_4 <- rps_out %>% filter(name %in% c("DDE_Tnp_1", "DDE_Tnp_4", "Plant_tran"))
MADF <- rps_out %>% filter(name %in% c("MADF_DNA_bdg", "MADF", "GT1", "SANT"))

potential_Harbinger <- rps_out %>% filter(qseqid %in% DDE_Tnp_1_4$aqeqid, qseqid %in% MADF$qseqid)

if(nrow(potential_Harbinger) > 0){
  
  potential_Harbinger <- rps_out %>%
    filter(qseiqd %in% potential_Harbinger, qlen < 8000) %>%
    dplyr::rename(seqnames = qseqid) %>%
    select(seqnames, qlen) %>%
    mutate(start = 1) %>%
    dplyr::rename(end = qlen) %>%
    base::unique() %>%
    as_granges()
  
  potential_Harbinger_seq <- getSeq(carp_seq, potential_Harbinger)
  names(potential_Harbinger_seq) <- seqnames(potential_Harbinger)
  writeXStringSet(potential_Harbinger_seq, paste0("CARP/", species_name, "/long_seq/potential_Harbinger.fa"))
  
}

# TcMariners
Tc1_DDE_HTH <- rps_out %>% filter(name %in% c("DDE_1", "DDE_3", "HTH_20", "HTH_23", "HTH_24", "HTH_28", "HTH_32",
                                              "HTH_5", "HTH_7", "HTH_AsnC-type", "HTH_Tnp_ISL3", "HTH_Tnp_Tc3_2"))

if(nrow(Tc1_DDE_HTH) > 0){
  
  potential_TcMar <- Tc1_DDE_HTH %>%
    dplyr::rename(seqnames = qseqid) %>%
    select(seqnames, qlen, strand) %>%
    mutate(start = 1) %>%
    dplyr::rename(end = qlen) %>%
    base::unique() %>%
    as_granges()
  
  potential_TcMar_seq <- getSeq(carp_seq, potential_TcMar)
  names(potential_TcMar_seq) <- seqnames(potential_TcMar)
  writeXStringSet(potential_TcMar_seq, paste0("CARP/", species_name, "/long_seq/potential_TcMar.fa"))
  
}

# hATs
Dimer_Tnp_hAT <- rps_out %>% filter(name == "Dimer_Tnp_hAT")

if(nrow(Dimer_Tnp_hAT) > 0){
  
  potential_hAT <- Dimer_Tnp_hAT %>%
    dplyr::rename(seqnames=qseqid) %>%
    select(seqnames, qlen, strand) %>%
    mutate(start = 1) %>%
    dplyr::rename(end = qlen) %>%
    base::unique() %>%
    as_granges()
  
  potential_hAT_seq <- getSeq(carp_seq, potential_hAT)
  names(potential_hAT_seq) <- seqnames(potential_hAT)
  writeXStringSet(potential_hAT_seq, paste0("CARP/", species_name, "/long_seq/potential_hAT.fa"))
  
}

# LTRs
LTRs <- rps_out %>%
  filter(grepl("RNase", name)) %>%
  group_by(qseqid) %>%
  dplyr::slice(1) %>%
  ungroup()


# Gypsys
Gypsy <- LTRs %>%
  filter(grepl("Ty3", name)) %>%
  dplyr::select(qseqid) %>%
  base::unique()

if(nrow(Gypsy) > 0){
  
  potential_Gypsys <- rps_out %>%
    filter(qseqid %in% Gypsy$qseqid) %>%
    dplyr::rename(seqnames = qseqid) %>%
    dplyr::select(seqnames, qlen, strand) %>%
    dplyr::rename(end = qlen) %>%
    dplyr::mutate(start = 1) %>%
    base::unique() %>%
    plyranges::as_granges()
  
  potential_Gypsys_seq <- getSeq(carp_seq, potential_Gypsys)
  names(potential_Gypsys_seq) <- seqnames(potential_Gypsys)
  writeXStringSet(potential_Gypsys_seq, paste0("CARP/", species_name, "/long_seq/potential_Gypsys.fa"))
  
}

# Copias
Copia <- LTRs %>%
  filter(grepl("Ty1", name)) %>%
  dplyr::select(qseqid) %>%
  base::unique()

if(nrow(Copia) > 0){
  potential_Copias <- rps_out %>%
    filter(qseqid %in% Copia$qseqid) %>%
    dplyr::rename(seqnames = qseqid) %>%
    dplyr::select(seqnames, qlen, strand) %>%
    dplyr::rename(end = qlen) %>%
    dplyr::mutate(start = 1) %>%
    base::unique() %>%
    plyranges::as_granges()
  
  potential_Copias_seq <- getSeq(carp_seq, potential_Copias)
  names(potential_Copias_seq) <- seqnames(potential_Copias)
  writeXStringSet(potential_Copias_seq, paste0("CARP/", species_name, "/long_seq/potential_Copias.fa"))
  
}

# DIRSs
DIRS <- LTRs %>%
  filter(grepl("DIRS", name)) %>%
  dplyr::select(qseqid) %>%
  base::unique()

if(nrow(DIRS) > 0){
  potential_DIRSs <- rps_out %>%
    filter(qseqid %in% DIRS$qseqid) %>%
    dplyr::rename(seqnames = qseqid) %>%
    dplyr::select(seqnames, qlen, strand) %>%
    dplyr::rename(end = qlen) %>%
    dplyr::mutate(start = 1) %>%
    base::unique() %>%
    plyranges::as_granges()
  
  potential_DIRSs_seq <- getSeq(carp_seq, potential_DIRSs)
  names(potential_DIRSs_seq) <- seqnames(potential_DIRSs)
  writeXStringSet(potential_DIRSs_seq, paste0("CARP/", species_name, "/long_seq/potential_DIRSs.fa"))
  
}

Bel <- LTRs %>%
  filter(grepl("Bel", name)) %>%
  dplyr::select(qseqid) %>%
  base::unique()

if(nrow(Bel) > 0){
  potential_Bels <- rps_out %>%
    filter(qseqid %in% Bel$qseqid) %>%
    dplyr::rename(seqnames = qseqid) %>%
    dplyr::select(seqnames, qlen, strand) %>%
    dplyr::rename(end = qlen) %>%
    dplyr::mutate(start = 1) %>%
    base::unique() %>%
    plyranges::as_granges()
  
  potential_Bels_seq <- getSeq(carp_seq, potential_Bels)
  names(potential_Bels_seq) <- seqnames(potential_Bels)
  writeXStringSet(potential_Bels_seq, paste0("CARP/", species_name, "/long_seq/potential_Bels.fa"))
  
}

blast_out <- read_tsv(system(paste0("blastn -task dc-megablast -num_threads 12 -query CARP/", species_name, "/long_seq/potential_Bels.fa -db ", genome_path, " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs")) %>%
  dplyr::filter(length >= 0.5 * qlen, pident >= 94) %>%
  dplyr::mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+"), # determine strand of hit
  start = case_when(strand == "+" ~ sstart, strand == "-" ~ send), # determine stranded start
  end = case_when(strand == "+" ~ send, strand == "-" ~ sstart))

to_mafft <- tibble(qseqid = names(table(blast_out$qseqid)), n = as.integer(table(blast_out$qseqid))) %>%
  filter(n >= 3)

repbase_blast_out <- read_tsv(system(paste0("blastn -task blastn -word_size 11 -num_threads 12 -query CARP/", species_name, "/long_seq/potential_Bels.fa -db ", repbase_path, " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs")) %>%
  group_by(qseqid) %>%
  dplyr::slice(1) %>%
  ungroup()

for(j in 1:nrow(to_mafft)){
  
  # subest base on query repeat
  to_align <- blast_out %>%
    filter(qseqid == to_mafft$qseqid[j]) %>%
    mutate(start = start - 2000, end = end + 2000,
           start = ifelse(start < 1, 1, start),
           end = ifelse(end > slen, slen, end)) %>%
    dplyr::slice(1:20) %>%
    as_granges()
  
  # get sequences from genome
  to_align_seq <- Biostrings::getSeq(genome_seq, to_align)
  
  # name sequences using scaffold_name:cordinates(strand) format
  names(to_align_seq) <- paste0(seqnames(to_align), ":", ranges(to_align), "(", strand(to_align), ")")
  
  # write compiled sequences to temporary fasta file
  Biostrings::writeXStringSet(x = to_align_seq, filepath = paste0("~/HT_Workflow/HT_Workflow/CARP/", species_name, "/long_seq/curation/temp.fa"), append = F)
  
  if(grepl(":", to_mafft$qseqid[j])){
    query_repeat_subbed <- sub(":.*", "", to_mafft$qseqid[j])
  } else {
    query_repeat_subbed <- sub("#.*", "", to_mafft$qseqid[j])
  }
  
  
  # create alignment of sequences, file named after original query with "/" character replaced with "_"
  system(paste0("mafft --localpair --thread 12 ~/HT_Workflow/HT_Workflow/CARP/", species_name, "/long_seq/curation/temp.fa > ~/HT_Workflow/HT_Workflow/CARP/", species_name, "/long_seq/curation/", species_name, "_", query_repeat_subbed, "_aligned.fa"))

}