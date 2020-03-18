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
recip_ranges <- recip_tibble %>%
  select(-hit_width, -names) %>%
  as_granges()

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

table(aipysurus_genome_annotation_tibble$type)

# create gffs and tibbles of each type
utr_gff <- aipysurus_genome_annotation_tibble %>%
  filter(type %in% c("three_prime_UTR", "five_prime_UTR")) %>%
  select(seqnames, start, end, strand, type, ID) %>%
  mutate(ID = sub(".utr.*", "", ID)) %>%
  inner_join(aipysurus_genome_annotation_seq_names) %>%
  select(-ID)
utr_gff_ranges <- as_granges(utr_gff)
mrna_gff <- aipysurus_genome_annotation_tibble %>%
  filter(type == "mRNA") %>%
  select(seqnames, start, end, strand, type, ID) %>%
  inner_join(aipysurus_genome_annotation_seq_names) %>%
  select(-ID)
mrna_gff_ranges <- as_granges(mrna_gff)
gene_gff <- aipysurus_genome_annotation_tibble %>%
  filter(type == "gene") %>%
  select(seqnames, start, end, strand, type, ID) %>%
  mutate(ID = sub("\\.TU\\.", ".model.", ID)) %>%
  inner_join(aipysurus_genome_annotation_seq_names) %>%
  select(-ID)
gene_gff_ranges <- as_granges(gene_gff)
cds_gff <- aipysurus_genome_annotation_tibble %>%
  filter(type == "CDS") %>%
  select(seqnames, start, end, strand, type, ID) %>%
  mutate(ID = sub("cds\\.", "", ID)) %>%
  inner_join(aipysurus_genome_annotation_seq_names) %>%
  select(-ID)
cds_gff_ranges <- as_granges(cds_gff)
exon_gff <- aipysurus_genome_annotation_tibble %>%
  filter(type == "exon") %>%
  mutate(ID = sub("\\.exon.*", "", ID)) %>%
  select(seqnames, start, end, strand, type, ID) %>%
  inner_join(aipysurus_genome_annotation_seq_names) %>%
  select(-ID)
exon_gff_ranges <- as_granges(exon_gff)

precede <- pair_precede(recip_ranges, aipysurus_genome_annotation) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames, -score, -phase, -X5_prime_partial, -X3_prime_partial, -Parent, -Name, -source) %>%
  rename(seqnames = granges.x.seqnames, repeat_start = granges.x.start, repeat_end = granges.x.end, ann_start = granges.y.start, ann_end = granges.y.end,
         repeat_strand = granges.x.strand, ann_strand = granges.y.strand) %>%
  mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand),
         ID = sub("\\.TU\\.", ".model.", ID), ID = sub("cds\\.", "", ID), ID = sub("\\.exon.*", "", ID), ID = sub(".utr.*", "", ID)) %>%
  inner_join(aipysurus_genome_annotation_seq_names) %>%
  select(-ID, -granges.y.width) %>%
  base::unique() %>%
  filter(grepl("UTR", type))


utr_distance_gappy <- pair_overlaps(recip_ranges, aipysurus_genome_annotation, maxgap = 1000) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames, -score, -phase, -X5_prime_partial, -X3_prime_partial, -Parent, -Name, -source) %>%
  rename(seqnames = granges.x.seqnames, repeat_start = granges.x.start, repeat_end = granges.x.end, ann_start = granges.y.start, ann_end = granges.y.end,
         repeat_strand = granges.x.strand, ann_strand = granges.y.strand) %>%
  mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand),
         ID = sub("\\.TU\\.", ".model.", ID), ID = sub("cds\\.", "", ID), ID = sub("\\.exon.*", "", ID), ID = sub(".utr.*", "", ID)) %>%
  inner_join(aipysurus_genome_annotation_seq_names) %>%
  select(-ID, -granges.y.width) %>%
  base::unique() %>%
  filter(grepl("UTR", type)) %>%
  # select(gene, type) %>%
  # base::unique() %>%
  View()

# covert ranges object to tibble, manipulate to remove unnecessary data and determine gene names
all_overlap_ranges_tibble <- all_overlap_ranges %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames, -score, -phase, -X5_prime_partial, -X3_prime_partial, -Parent, -Name, -source) %>%
  rename(seqnames = granges.x.seqnames, repeat_start = granges.x.start, repeat_end = granges.x.end, ann_start = granges.y.start, ann_end = granges.y.end,
         repeat_strand = granges.x.strand, ann_strand = granges.y.strand) %>%
  mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand),
         ID = sub("\\.TU\\.", ".model.", ID), ID = sub("cds\\.", "", ID), ID = sub("\\.exon.*", "", ID), ID = sub(".utr.*", "", ID)) %>%
  inner_join(aipysurus_genome_annotation_seq_names) %>%
  select(-ID, -granges.y.width, -granges.x.width) %>%
  base::unique()

# write tibble to file
write_tsv(all_overlap_ranges_tibble, "GeneInteraction/insertion_sites.tsv", col_names = T)

# extract repeat insertions within 5000bp of 5' UTRs
utr_distance_gappy <- pair_overlaps(recip_ranges, aipysurus_genome_annotation, maxgap = 5000) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames, -score, -phase, -X5_prime_partial, -X3_prime_partial, -Parent, -Name, -source) %>%
  rename(seqnames = granges.x.seqnames, repeat_start = granges.x.start, repeat_end = granges.x.end, ann_start = granges.y.start, ann_end = granges.y.end,
         repeat_strand = granges.x.strand, ann_strand = granges.y.strand) %>%
  mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand),
         ID = sub("\\.TU\\.", ".model.", ID), ID = sub("cds\\.", "", ID), ID = sub("\\.exon.*", "", ID), ID = sub(".utr.*", "", ID)) %>%
  inner_join(aipysurus_genome_annotation_seq_names) %>%
  select(-ID, -granges.y.width, -granges.x.width) %>%
  base::unique() %>%
  filter(grepl("UTR", type))

# filter repeat insertions up to be upstream of 5' UTRs
upstream <- utr_distance_gappy %>%
  filter(type == "five_prime_UTR") %>%
  filter((ann_strand == "+" & repeat_end < ann_start) | (ann_strand == "-" & repeat_start > ann_end))
# write_tsv(upstream, "GeneInteraction/upstream.tsv", col_names = T)

# select insertions into mRNA
mrna_overlaps <- all_overlap_ranges_tibble %>%
  filter(type == "mRNA")
write_tsv(mrna_overlaps, "GeneInteraction/mrna_overlaps.tsv", col_names = T)
