library(BSgenome)
library(GenomicRanges)
library(plyranges)
library(tidyverse)

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
  select(seqnames, start, end, strand, ID) %>%
  mutate(ID = sub(".utr.*", "", ID)) %>%
  inner_join(aipysurus_genome_annotation_seq_names) %>%
  select(-ID) %>%
  base::unique()
utr_gff_ranges <- as_granges(utr_gff)

mrna_gff <- aipysurus_genome_annotation_tibble %>%
  filter(type == "mRNA") %>%
  select(seqnames, start, end, strand, ID) %>%
  inner_join(aipysurus_genome_annotation_seq_names) %>%
  select(-ID) %>%
  base::unique()
mrna_gff_ranges <- as_granges(mrna_gff)

gene_gff <- aipysurus_genome_annotation_tibble %>%
  filter(type == "gene") %>%
  select(seqnames, start, end, strand, ID) %>%
  mutate(ID = sub("\\.TU\\.", ".model.", ID)) %>%
  inner_join(aipysurus_genome_annotation_seq_names) %>%
  select(-ID) %>%
  base::unique()
gene_gff_ranges <- as_granges(gene_gff)

cds_gff <- aipysurus_genome_annotation_tibble %>%
  filter(type == "CDS") %>%
  select(seqnames, start, end, strand, ID) %>%
  mutate(ID = sub("cds\\.", "", ID)) %>%
  inner_join(aipysurus_genome_annotation_seq_names) %>%
  select(-ID) %>%
  base::unique()
cds_gff_ranges <- as_granges(cds_gff)

exon_gff <- aipysurus_genome_annotation_tibble %>%
  filter(type == "exon") %>%
  mutate(ID = sub("\\.exon.*", "", ID)) %>%
  select(seqnames, start, end, strand, ID) %>%
  inner_join(aipysurus_genome_annotation_seq_names) %>%
  select(-ID) %>%
  base::unique()
exon_gff_ranges <- as_granges(exon_gff)

# Read in repeat annotation
aipysurus_repeat_annotation <- read_table("GeneInteraction/kmer_49.pilon_x2.sorted.fasta.out", skip = 3,
                                          col_names = c("sw_score", "perc_div", "perc_del", "perc_ins", "seqnames", "start", "end", "remain", "strand",
                                                        "m_repeat", "repeat_type", "r_begin", "r_end", "r_left", "repeat_no", "overlap")) 

aipysurus_repeat_annotation_corrected <- aipysurus_repeat_annotation %>%
  mutate(strand = sub("C", "-", strand)) %>%
  filter(start < end, strand %in% c("-", "+"), overlap != "*", !repeat_type %in% c("Unknown", "Simple_repeat")) %>%
  mutate(full_repeat = paste0(m_repeat, "#", repeat_type)) %>%
  dplyr::select(seqnames, start, end, strand, full_repeat, perc_div)

aipysurus_lines_annotation <- aipysurus_repeat_annotation_corrected %>%
  filter(grepl("#LINE", full_repeat), end - start > 50, perc_div <= 25)


# extract LINEs, limit by size and divergence
aipysurus_lines_annotation <- aipysurus_repeat_annotation_corrected %>%
  filter(grepl("#LINE", full_repeat), !grepl("Penelope", full_repeat), end - start > 50, perc_div <= 25)
aipysurus_lines_annotation_ranges <- plyranges::as_granges(aipysurus_lines_annotation)

# find overlaps of genes and lines
genes_lines_overlap <- pair_overlaps(gene_gff_ranges, aipysurus_lines_annotation_ranges) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames) %>%
  rename(seqnames = granges.x.seqnames, ann_start = granges.x.start, ann_end = granges.x.end, repeat_start = granges.y.start, repeat_end = granges.y.end,
         ann_strand = granges.x.strand, repeat_strand = granges.y.strand) %>%
  mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand))

# extract genes containing LINEs
genes_contain_lines <- genes_lines_overlap %>%
  select(seqnames, ann_start, ann_end, ann_strand, gene) %>%
  base::unique()
genes_contain_lines_ranges <- genes_contain_lines %>%
  dplyr::rename(start = ann_start, end = ann_end, strand = ann_strand) %>%
  plyranges::as_granges()

# extract LINEs within genes
lines_in_genes <- genes_lines_overlap %>%
  dplyr::select(seqnames, repeat_start, repeat_end, repeat_strand, full_repeat, perc_div) %>%
  base::unique()
lines_in_genes_ranges <- lines_in_genes %>%
  dplyr::rename(start = repeat_start, end = repeat_end, strand = repeat_strand) %>%
  plyranges::as_granges()

# extract LINEs in exons
lines_exons_overlap <- pair_overlaps(exon_gff_ranges, lines_in_genes_ranges) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames) %>%
  rename(seqnames = granges.x.seqnames, ann_start = granges.x.start, ann_end = granges.x.end, repeat_start = granges.y.start, repeat_end = granges.y.end,
         ann_strand = granges.x.strand, repeat_strand = granges.y.strand) %>%
  mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand))

lines_exons_overlap_ranges <- lines_exons_overlap %>%
  rename(start = repeat_start, end = repeat_end, strand = repeat_strand) %>%
  filter(!grepl("LOC", gene), !grepl("\\.", gene)) %>%
  select(seqnames, start, end, strand, full_repeat, gene) %>%
  plyranges::as_granges() %>%
  plyranges::reduce_ranges_directed(full_repeat = unique(full_repeat), gene = unique(gene))
  
lines_exons_overlap_refined <- lines_exons_overlap_ranges %>%
  as_tibble() %>%
  mutate(full_repeat = as.character(full_repeat), gene = as.character(gene), seqnames = as.character(seqnames), strand = as.character(strand))

lines_exons_overlap_ranges <- plyranges::as_granges(lines_exons_overlap_refined)

# extract LINEs not in genes
lines_nongene_ranges <- filter_by_non_overlaps(aipysurus_lines_annotation_ranges, lines_in_genes_ranges)

lines_nongene_ranges %>% GenomicRanges::reduce(min.gapwidth = 500)

lines_nongene_ranges %>% 
  reduce_ranges_directed(full_repeat = base::unique(full_repeat), perc_div = base::unique(perc_div)) %>%
  dplyr::mutate(perc_div = base::unlist(perc_div)) %>%
  dplyr::mutate(full_repeat = base::unlist(full_repeat))

# extract LINEs in UTRs
lines_utr_overlap <- pair_overlaps(utr_gff_ranges, lines_in_genes_ranges) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames) %>%
  filter(!grepl("LOC", gene), !grepl("\\.", gene)) %>%
  dplyr::rename(seqnames = granges.x.seqnames, ann_start = granges.x.start, ann_end = granges.x.end, repeat_start = granges.y.start, repeat_end = granges.y.end,
         ann_strand = granges.x.strand, repeat_strand = granges.y.strand) %>%
  dplyr::mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand))

lines_in_utr_ranges <- lines_utr_overlap %>%
  dplyr::rename(start = repeat_start, end = repeat_end, strand = repeat_strand) %>%
  dplyr::select(seqnames, start, end, strand, gene, full_repeat, perc_div) %>%
  plyranges::as_granges() 


lines_in_utr_reduced_ranges <- lines_in_utr_ranges %>%
  reduce_ranges_directed(gene = base::unique(gene), full_repeat = base::unique(full_repeat), perc_div = base::unique(perc_div)) %>%
  dplyr::mutate(perc_div = base::unlist(perc_div)) %>%
  dplyr::mutate(full_repeat = base::unlist(full_repeat)) %>%
  dplyr::mutate(gene = base::unlist(gene))

lines_in_utr_seq <- Biostrings::getSeq(genome_seq, lines_in_utr_reduced_ranges)
names(lines_in_utr_seq) <- paste0(lines_in_utr_reduced_ranges$gene, "#", lines_in_utr_reduced_ranges$full_repeat)
writeXStringSet(lines_in_utr_seq, "GeneInteraction/utr_repeats.fasta")

lines_in_coding_exons_ranges <- filter_by_non_overlaps(lines_exons_overlap_ranges, lines_in_utr_ranges)

lines_in_coding_exons_seq <- Biostrings::getSeq(genome_seq, lines_in_coding_exons_ranges)
names(lines_in_coding_exons_seq) <- paste0(lines_in_coding_exons_ranges$gene, "#", lines_in_coding_exons_ranges$full_repeat)
writeXStringSet(lines_in_coding_exons_seq, "GeneInteraction/coding_exon_repeats.fasta")


# # if going to reduce
# r <- GenomicRanges::reduce(lines_nongene_ranges, with.revmap = TRUE, min.gapwidth = 500)
# revmap_val <- mcols(r)$revmap
# mcols(r) <- DataFrame(full_repeat = base::unique(extractList(mcols(lines_nongene_ranges)$full_repeat, revmap_val)))
# 
# lines_nongene_reduced <- as_tibble(r) %>%
#   mutate(full_repeat = as.character(full_repeat))
# 
# lines_nongene_reduced %>%
#   filter(grepl("Snek", full_repeat)) %>%
#   View()
