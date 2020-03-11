# Script for XXX

# setwd("~/Analysis/Snake/HT_Workflow/")

library(BSgenome)
library(GenomicRanges)
library(plyranges)
library(tidyverse)


# Set repeat and species names
repeat_list <- tibble(repeat_name = c("Proto2-Snek", "Rex1-Snek_1", "Rex1-Snek_2", "Rex1-Snek_3", "Rex1-Snek_4", "RTE-Snek"))

# set variables
genome_dir <- "~/Genomes/Reptiles/"
species_name <- "Aipysurus_laevis"
genome_name <- "kmer_49.pilon_x2.sorted.fasta"
print(species_name)

# set and read in genome and index
genome_path <- paste0(genome_dir, "/", species_name, "/", genome_name)
genome_seq <- Biostrings::readDNAStringSet(filepath = genome_path)
gc()
names(genome_seq) <- sub(" .*", "", names(genome_seq))
genome_idx <- tibble(seqnames = names(genome_seq), width = as.double(width(genome_seq)))

# Search for lines in genome
blast_out <- read_tsv(system(paste0("blastn -dust yes -query compiled_lines.fasta -db ", genome_path, " -outfmt \"6 sseqid sstart send pident qcovs bitscore length mismatch evalue qseqid\""), intern = TRUE), col_names = c("seqnames", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "qseqid")) %>%
  filter(length > 50, pident >= 97) %>%
  mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+"),
         start = case_when(sstart < send ~ sstart, send < sstart ~ send),
         end = case_when(sstart > send ~ sstart, send > sstart ~ send))

# Create bed of blast output
bed <- blast_out %>%
  inner_join(genome_idx) %>%
  dplyr::select(seqnames, start, end, strand) %>%
  as_granges() %>%
  GenomicRanges::reduce(min.gapwidth = 500)
bed$names <- paste0(seqnames(bed), ":", ranges(bed), "(", strand(bed), ")")

# Get seqs of blast hits and write to file
seqs <- Biostrings::getSeq(genome_seq, bed)
names(seqs) <- bed$names
Biostrings::writeXStringSet(x = seqs, filepath = "GeneInteraction/initialHits.fa")

# align hits to original search and select best hit for each based on bitscore
blastn_recip <- read_tsv(system(paste0("blastn -dust yes -query GeneInteraction/initialHits.fa -subject compiled_lines.fasta  -outfmt \"6 qseqid sseqid qstart qend sstart send pident qcovs bitscore length mismatch evalue\""), intern = TRUE), col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue")) %>%
  group_by(qseqid) %>%
  dplyr::slice(1) %>%
  select(qseqid, sseqid) %>%
  dplyr::rename(names = qseqid, repeat_name = sseqid) %>%
  ungroup()

recip_tibble <- bed %>%
  as_tibble() %>%
  mutate(seqnames = as.character(seqnames)) %>%
  dplyr::rename(hit_width = width) %>%
  inner_join(blastn_recip)

recip_ranges <- as_granges(recip_tibble)

# read in aipysurus annotation gff and fasta
aipysurus_genome_annotation <- plyranges::read_gff(file = "~/Genomes/Reptiles/Aipysurus_laevis/kmer_49.pilon_x2.sorted_annotation.gff")
aipysurus_genome_annotation_seq <- Biostrings::readDNAStringSet(filepath = "~/Genomes/Reptiles/Aipysurus_laevis/kmer_49.pilon_x2.sorted_annotation.fna")
gc()

# Split annotation faster header into ID and name
aipysurus_genome_annotation_seq_names <- tibble(names = names(aipysurus_genome_annotation_seq)) %>%
  mutate(ID = sub("_([^x]*)$", "", names), gene = sub(".*_", "", names)) %>%
  select(-names)

# convert annotation gff to tibble
aipysurus_genome_annotation_tibble <- aipysurus_genome_annotation %>%
  as_tibble()

# create gffs and tibbles of each type
utr_gff <- aipysurus_genome_annotation_tibble %>%
  filter(type %in% c("three_prime_UTR", "five_prime_UTR")) %>%
  select(seqnames, start, end, strand, type, ID) %>%
  mutate(ID = sub(".utr.*", "", ID)) %>%
  inner_join(aipysurus_genome_annotation_seq_names)
utr_gff_ranges <- as_granges(utr_gff)
mrna_gff <- aipysurus_genome_annotation_tibble %>%
  filter(type == "mRNA") %>%
  select(seqnames, start, end, strand, type, ID) %>%
  inner_join(aipysurus_genome_annotation_seq_names)
mrna_gff_ranges <- as_granges(mrna_gff)
gene_gff <- aipysurus_genome_annotation_tibble %>%
  filter(type == "gene") %>%
  select(seqnames, start, end, strand, type, ID) %>%
  mutate(ID = sub("\\.TU\\.", ".model.", ID)) %>%
  inner_join(aipysurus_genome_annotation_seq_names)
gene_gff_ranges <- as_granges(gene_gff)
cds_gff <- aipysurus_genome_annotation_tibble %>%
  filter(type == "CDS") %>%
  select(seqnames, start, end, strand, type, ID) %>%
  mutate(ID = sub("cds\\.", "", ID)) %>%
  inner_join(aipysurus_genome_annotation_seq_names)
cds_gff_ranges <- as_granges(cds_gff)
exon_gff <- aipysurus_genome_annotation_tibble %>%
  filter(type == "exon") %>%
  mutate(ID = sub("\\.exon.*", "", ID)) %>%
  select(seqnames, start, end, strand, type, ID) %>%
  inner_join(aipysurus_genome_annotation_seq_names)
exon_gff_ranges <- as_granges(exon_gff)

# intersect ranges
utr_overlap_ranges <- suppressWarnings(find_overlaps(recip_ranges, utr_gff_ranges))
mrna_overlap_ranges <- suppressWarnings(find_overlaps(recip_ranges, mrna_gff_ranges))
cds_overlap_ranges <- suppressWarnings(find_overlaps(recip_ranges, cds_gff_ranges))
exon_overlap_ranges <- suppressWarnings(find_overlaps(recip_ranges, exon_gff_ranges))

