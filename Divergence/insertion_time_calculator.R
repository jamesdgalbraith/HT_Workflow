# Script for plotting Jukes Cantor divergence of repeats found in sea snakes

library(BSgenome)
library(GenomicRanges)
library(plyranges)
library(tidyverse)

# create scratch directory
if(!dir.exists(paste0("~/HT_Workflow/HT_Workflow/Divergence/scratch/"))){
  dir.create(paste0("~/HT_Workflow/HT_Workflow/Divergence/scratch/"), recursive = T)
}

# make list of repeat names
repeat_list <- tibble(repeat_name = c("Proto2-Snek", "Rex1-Snek_1", "Rex1-Snek_2", "Rex1-Snek_3", "Rex1-Snek_4", "RTE-Snek_1", "RTE-Snek_2"))

# list species
species_list <- tibble(species = c("Aipysurus_laevis", "Hydrophis_curtis", "Emydocephalus_ijimae"),
                       genome = c("kmer_49.pilon_x2.sorted.fasta", "Hcur1.v1.1.fasta", "emyIji_1.0.fasta"))

i=1

# set genome path
genome_path <- paste0("~/Genomes/Reptiles/", species_list$species[i], "/", species_list$genome[i])

# read in genome and garbage collect
genome_seq <- Biostrings::readDNAStringSet(filepath = genome_path)
gc()
names(genome_seq) <- sub(" .*", "", names(genome_seq))

# search for repeats, require hits to be >100bp and 90% pident or higher
blastn <- read_tsv(system(paste0("blastn -dust yes -query ~/HT_Workflow/HT_Workflow/compiled_LINEs.fasta -db ", genome_path, " -outfmt \"6 sseqid sstart send pident qcovs bitscore length mismatch evalue qseqid\""), intern = TRUE), col_names = c("seqnames", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "qseqid")) %>%
  as_tibble() %>%
  dplyr::filter(length >= 100, pident >= 90) %>%
  mutate(pident_rnd = ceiling(pident), div = 100 - pident, d = mismatch/length, seqnames = as.character(seqnames), qseqid = as.character(qseqid),
         jc_dist = (-3 / 4) * log(1 - (4 * d / 3))) # calculate Jukes-Cantor distance

# convert output to ranges and reduce
blastn_ranges <- blastn %>%
  dplyr::mutate(strand = ifelse(sstart < send, "+", "-")) %>%
  dplyr::mutate(start = ifelse(strand == "+", sstart, send),
                end = ifelse(strand == "-", sstart, send)) %>%
  dplyr::select(seqnames, start, end, strand) %>%
  plyranges::as_granges() %>%
  GenomicRanges::reduce()

# get seq and write to file
blastn_seq <- BSgenome::getSeq(genome_seq, blastn_ranges)
names(blastn_seq) <- paste0(seqnames(blastn_ranges), ":", ranges(blastn_ranges), "(", strand(blastn_ranges), ")")
writeXStringSet(blastn_seq, paste0("~/HT_Workflow/HT_Workflow/Divergence/scratch/", species_list$species[i], "blastn_seq.fa"))

# search hits back against repeat library
recip_blastn <- read_tsv(system(paste0("blastn -dust yes -query ~/HT_Workflow/HT_Workflow/Divergence/scratch/", species_list$species[i], "blastn_seq.fa -subject ~/HT_Workflow/HT_Workflow/compiled_LINEs.fasta -outfmt \"6 qseqid sseqid sstart send pident qcovs bitscore length mismatch evalue\""), intern = TRUE), col_names = c("qseqid", "seqnames", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue")) %>%
  as_tibble() %>%
  dplyr::filter(length >= 100)

# select top hits and calculate insertion time based on substitution per 10 year generation (in Ludington and Sanders, in review at Molecular Ecology)
top_hits <- recip_blastn %>%
  group_by(qseqid) %>%
  arrange(-pident) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  mutate(sps = mismatch/length, gen = sps/1.25e-08, time = gen *10)

# calculate divergence time
div_time <- top_hits %>%
  group_by(seqnames) %>%
  summarise(mya = mean(time)/1e6, n = n(), sd = sd(time)/1e6, var = var(time))

# write tsv
write_tsv(div_time, paste0("~/HT_Workflow/HT_Workflow/Divergence/scratch/", species_list$species[i], "_div_time.tsv"))

# create densoty plot of pident
ggplot(top_hits, aes(x = pident, colour = seqnames)) + geom_density(alpha = 0.1)
ggplot(top_hits, aes(x = length, colour = seqnames)) + geom_density(alpha = 0.1)
