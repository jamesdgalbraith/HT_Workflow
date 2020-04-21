library(BSgenome)
library(GenomicRanges)
library(plyranges)
library(tidyverse)


# set and read in genome and index
genome_dir <- "~/Genomes/Reptiles/"
species_name <- "Aipysurus_laevis"
genome_name <- "kmer_49.pilon_x2.sorted.fasta"

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

# Read in repeat annotation
aipysurus_repeat_annotation <- readr::read_table("GeneInteraction/kmer_49.pilon_x2.sorted.fasta.out", skip = 3,
                                          col_names = c("sw_score", "perc_div", "perc_del", "perc_ins", "seqnames", "start", "end", "remain", "strand",
                                                        "m_repeat", "repeat_type", "r_begin", "r_end", "r_left", "repeat_no", "overlap")) 

aipysurus_repeat_annotation_corrected <- aipysurus_repeat_annotation %>%
  dplyr::mutate(strand = sub("C", "-", strand)) %>%
  dplyr::filter(start < end, strand %in% c("-", "+"), !repeat_type %in% c("Unknown", "Simple_repeat"), perc_div <= 25) %>%
  # mutate(full_repeat = paste0(m_repeat, "#", repeat_type)) %>%
  dplyr::select(seqnames, start, end, strand, m_repeat, repeat_type, perc_div, overlap)

# create gffs and tibbles of each type
utr_gff <- aipysurus_genome_annotation_tibble %>%
  dplyr::filter(type %in% c("three_prime_UTR", "five_prime_UTR")) %>%
  dplyr::select(seqnames, start, end, strand, ID, type) %>%
  dplyr::mutate(ID = sub(".utr.*", "", ID), type = as.character(type)) %>%
  dplyr::inner_join(aipysurus_genome_annotation_seq_names) %>%
  dplyr::select(-ID) %>%
  base::unique()

utr_gff_ranges <- plyranges::as_granges(utr_gff)

five_prime_utr_gff_ranges <- utr_gff_ranges %>%
  dplyr::filter(type == "five_prime_UTR")

mrna_gff <- aipysurus_genome_annotation_tibble %>%
  dplyr::filter(type == "mRNA") %>%
  dplyr::select(seqnames, start, end, strand, ID) %>%
  dplyr::inner_join(aipysurus_genome_annotation_seq_names) %>%
  dplyr::select(-ID) %>%
  base::unique()
mrna_gff_ranges <- plyranges::as_granges(mrna_gff)

gene_gff <- aipysurus_genome_annotation_tibble %>%
  dplyr::filter(type == "gene") %>%
  dplyr::select(seqnames, start, end, strand, ID) %>%
  dplyr::mutate(ID = sub("\\.TU\\.", ".model.", ID)) %>%
  dplyr::inner_join(aipysurus_genome_annotation_seq_names) %>%
  dplyr::select(-ID) %>%
  base::unique()
gene_gff_ranges <- plyranges::as_granges(gene_gff)

cds_gff <- aipysurus_genome_annotation_tibble %>%
  dplyr::filter(type == "CDS") %>%
  dplyr::select(seqnames, start, end, strand, ID) %>%
  dplyr::mutate(ID = sub("cds\\.", "", ID)) %>%
  dplyr::inner_join(aipysurus_genome_annotation_seq_names) %>%
  dplyr::select(-ID) %>%
  base::unique()
cds_gff_ranges <- plyranges::as_granges(cds_gff)

exon_gff <- aipysurus_genome_annotation_tibble %>%
  dplyr::filter(type == "exon") %>%
  dplyr::mutate(ID = sub("\\.exon.*", "", ID)) %>%
  dplyr::select(seqnames, start, end, strand, ID) %>%
  dplyr::inner_join(aipysurus_genome_annotation_seq_names) %>%
  dplyr::select(-ID) %>%
  base::unique()
exon_gff_ranges <- plyranges::as_granges(exon_gff)

# extract LINEs, limit by size and divergence
aipysurus_lines_annotation_ranges <- aipysurus_repeat_annotation_corrected %>%
  dplyr::filter(grepl("LINE", repeat_type), !grepl("Penelope", repeat_type), end - start > 50, overlap != "*") %>%
  plyranges::as_granges() 
aipysurus_lines_annotation_ranges_r <- GenomicRanges::reduce(aipysurus_lines_annotation_ranges, with.revmap = TRUE, min.gapwidth = 250)
revmap_val <- mcols(aipysurus_lines_annotation_ranges_r)$revmap
mcols(aipysurus_lines_annotation_ranges_r) <- DataFrame(repeat_type = extractList(mcols(aipysurus_lines_annotation_ranges)$repeat_type, revmap_val))
aipysurus_lines_annotation_ranges_r <- aipysurus_lines_annotation_ranges_r %>%
  dplyr::mutate(repeat_type = base::unique(unlist(repeat_type)))

aipysurus_lines_annotation_ranges <- aipysurus_lines_annotation_ranges_r %>%
  dplyr::filter(seqnames %in% seqnames(gene_gff_ranges) | seqnames %in% seqnames(exon_gff_ranges) | seqnames %in% seqnames(five_prime_utr_gff_ranges))

aipysurus_sneks_annotation_ranges <- aipysurus_repeat_annotation_corrected %>%
  dplyr::filter(grepl("Snek", m_repeat), end - start > 50) %>%
  dplyr::select(seqnames, start, end, strand, m_repeat) %>%
  plyranges::as_granges()

aipysurus_sneks_annotation_ranges <- aipysurus_sneks_annotation_ranges %>%
  dplyr::filter(seqnames %in% seqnames(gene_gff_ranges) | seqnames %in% seqnames(exon_gff_ranges) | seqnames %in% seqnames(five_prime_utr_gff_ranges))

# find overlaps of genes and lines
genes_lines_overlap <- plyranges::pair_overlaps(gene_gff_ranges, aipysurus_lines_annotation_ranges) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames) %>%
  dplyr::rename(seqnames = granges.x.seqnames, ann_start = granges.x.start, ann_end = granges.x.end, repeat_start = granges.y.start, repeat_end = granges.y.end,
         ann_strand = granges.x.strand, repeat_strand = granges.y.strand) %>%
  dplyr::mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand))

genes_sneks_overlap <- plyranges::pair_overlaps(gene_gff_ranges, aipysurus_sneks_annotation_ranges) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames) %>%
  dplyr::rename(seqnames = granges.x.seqnames, ann_start = granges.x.start, ann_end = granges.x.end, repeat_start = granges.y.start, repeat_end = granges.y.end,
         ann_strand = granges.x.strand, repeat_strand = granges.y.strand) %>%
  dplyr::mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand))

# extract genes containing LINEs
genes_contain_lines <- genes_lines_overlap %>%
  dplyr::select(seqnames, ann_start, ann_end, ann_strand, gene) %>%
  base::unique()
genes_contain_lines_ranges <- genes_contain_lines %>%
  dplyr::rename(start = ann_start, end = ann_end, strand = ann_strand) %>%
  plyranges::as_granges()

genes_contain_sneks <- genes_sneks_overlap %>%
  dplyr::select(seqnames, ann_start, ann_end, ann_strand, gene, m_repeat) %>%
  base::unique()
genes_contain_sneks_ranges <- genes_contain_sneks %>%
  dplyr::rename(start = ann_start, end = ann_end, strand = ann_strand) %>%
  plyranges::as_granges()

# extract LINEs within genes
lines_in_genes <- genes_lines_overlap %>%
  dplyr::select(seqnames, repeat_start, repeat_end, repeat_strand, repeat_type, gene) %>%
  base::unique()

lines_in_genes_ranges <- lines_in_genes %>%
  dplyr::rename(start = repeat_start, end = repeat_end, strand = repeat_strand) %>%
  plyranges::as_granges()

sneks_in_genes <- genes_sneks_overlap %>%
  dplyr::select(seqnames, repeat_start, repeat_end, repeat_strand, m_repeat, gene) %>%
  base::unique()

sneks_in_genes_ranges <- sneks_in_genes %>%
  dplyr::rename(start = repeat_start, end = repeat_end, strand = repeat_strand) %>%
  plyranges::as_granges()

# extract LINEs in exons
lines_exons_overlap <- pair_overlaps(exon_gff_ranges, aipysurus_lines_annotation_ranges) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames) %>%
  dplyr::rename(seqnames = granges.x.seqnames, ann_start = granges.x.start, ann_end = granges.x.end, repeat_start = granges.y.start, repeat_end = granges.y.end,
         ann_strand = granges.x.strand, repeat_strand = granges.y.strand) %>%
  dplyr::mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand)) %>%
  base::unique()

lines_in_exons <- lines_exons_overlap %>%
  dplyr::select(seqnames, repeat_start, repeat_end, repeat_strand, repeat_type, gene)

lines_in_exons_ranges <- lines_in_exons %>%
  dplyr::rename(start = repeat_start, end = repeat_end, strand = repeat_strand) %>%
  dplyr::select(seqnames, start, end, strand, repeat_type, gene) %>%
  plyranges::as_granges()
  
sneks_exons_overlap <- pair_overlaps(exon_gff_ranges, aipysurus_sneks_annotation_ranges) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames) %>%
  dplyr::rename(seqnames = granges.x.seqnames, ann_start = granges.x.start, ann_end = granges.x.end, repeat_start = granges.y.start, repeat_end = granges.y.end,
         ann_strand = granges.x.strand, repeat_strand = granges.y.strand) %>%
  dplyr::mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand)) %>%
  base::unique()

sneks_in_exons <- sneks_exons_overlap %>%
  dplyr::select(seqnames, repeat_start, repeat_end, repeat_strand, m_repeat, gene)

sneks_in_exons_ranges <- sneks_in_exons %>%
  dplyr::rename(start = repeat_start, end = repeat_end, strand = repeat_strand) %>%
  dplyr::select(seqnames, start, end, strand, m_repeat, gene) %>%
  plyranges::as_granges()

# extract LINEs in UTRs
lines_utr_overlap <- pair_overlaps(utr_gff_ranges, lines_in_genes_ranges) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames) %>%
  dplyr::rename(seqnames = granges.x.seqnames, ann_start = granges.x.start, ann_end = granges.x.end, repeat_start = granges.y.start, repeat_end = granges.y.end,
         ann_strand = granges.x.strand, repeat_strand = granges.y.strand, gene = gene.x) %>%
  dplyr::mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand)) %>%
  base::unique()

lines_in_utr_ranges <- lines_utr_overlap %>%
  dplyr::rename(start = repeat_start, end = repeat_end, strand = repeat_strand) %>%
  dplyr::select(seqnames, start, end, strand, gene, repeat_type) %>%
  plyranges::as_granges() 

sneks_utr_overlap <- pair_overlaps(utr_gff_ranges, sneks_in_genes_ranges) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames) %>%
  dplyr::rename(seqnames = granges.x.seqnames, ann_start = granges.x.start, ann_end = granges.x.end, repeat_start = granges.y.start, repeat_end = granges.y.end,
                ann_strand = granges.x.strand, repeat_strand = granges.y.strand, gene = gene.x) %>%
  dplyr::mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand)) %>%
  base::unique()

sneks_in_utr_ranges <- sneks_utr_overlap %>%
  dplyr::rename(start = repeat_start, end = repeat_end, strand = repeat_strand) %>%
  dplyr::select(seqnames, start, end, strand, gene, m_repeat) %>%
  plyranges::as_granges() 

# setract lines in coding regions
lines_in_coding_exons_ranges <- filter_by_non_overlaps(lines_in_exons_ranges, lines_in_utr_ranges) %>%
  as_tibble() %>%
  base::unique() %>%
  plyranges::as_granges()

sneks_in_coding_exons_ranges <- filter_by_non_overlaps(sneks_in_exons_ranges, sneks_in_utr_ranges) %>%
  as_tibble() %>%
  base::unique() %>%
  plyranges::as_granges()

# extract LINEs not in genes
lines_near_genes_ranges <- aipysurus_lines_annotation_ranges %>%
  dplyr::filter(seqnames %in% seqnames(gene_gff_ranges)) %>%
  base::unique()

lines_nongene_ranges <- filter_by_non_overlaps(lines_near_genes_ranges, gene_gff_ranges) %>%
  dplyr::filter(seqnames %in% aipysurus_genome_annotation_tibble$seqnames)

sneks_near_genes_ranges <- aipysurus_sneks_annotation_ranges %>%
  dplyr::filter(seqnames %in% seqnames(gene_gff_ranges)) %>%
  base::unique()

sneks_nongene_ranges <- filter_by_non_overlaps(sneks_near_genes_ranges, gene_gff_ranges) %>%
  dplyr::filter(seqnames %in% aipysurus_genome_annotation_tibble$seqnames)

# find LINEs within 5000bp of utrs
utr_5000_to_lines <- pair_overlaps(lines_nongene_ranges, five_prime_utr_gff_ranges, maxgap = 5000) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames) %>%
  dplyr::rename(seqnames = granges.x.seqnames, ann_start = granges.y.start, ann_end = granges.y.end, repeat_start = granges.x.start, repeat_end = granges.x.end,
         ann_strand = granges.y.strand, repeat_strand = granges.x.strand) %>%
  dplyr::mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand)) %>%
  base::unique() %>%
  dplyr::filter((ann_strand == "-" & repeat_start > ann_end) | (ann_strand == "+" & repeat_end < ann_start))

utr_5000_to_sneks <- pair_overlaps(sneks_nongene_ranges, five_prime_utr_gff_ranges, maxgap = 5000) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames) %>%
  dplyr::rename(seqnames = granges.x.seqnames, ann_start = granges.y.start, ann_end = granges.y.end, repeat_start = granges.x.start, repeat_end = granges.x.end,
         ann_strand = granges.y.strand, repeat_strand = granges.x.strand) %>%
  dplyr::mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand)) %>%
  base::unique() %>%
  dplyr::filter((ann_strand == "-" & repeat_start > ann_end) | (ann_strand == "+" & repeat_end < ann_start))

# select sneks upstream of 5' utrs
sneks_upstream <- pair_overlaps(aipysurus_sneks_annotation_ranges, five_prime_utr_gff_ranges, maxgap = 5000) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames) %>%
  dplyr::rename(seqnames = granges.x.seqnames, ann_start = granges.y.start, ann_end = granges.y.end, repeat_start = granges.x.start, repeat_end = granges.x.end,
         ann_strand = granges.y.strand, repeat_strand = granges.x.strand) %>%
  dplyr::mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand)) %>%
  base::unique() %>%
  dplyr::filter((ann_strand == "-" & repeat_start > ann_end) | (ann_strand == "+" & repeat_end < ann_start)) %>%
  dplyr::select(seqnames, repeat_start, repeat_end, repeat_strand, m_repeat, gene) %>%
  dplyr::rename(start = repeat_start, end = repeat_end, strand = repeat_strand) %>%
  plyranges::as_granges()

# remove upstream sneks found in genes
sneks_upstream_nongene <- filter_by_non_overlaps(sneks_upstream, gene_gff_ranges) %>%
  as_tibble() %>%
  dplyr::arrange(gene) %>%
  select(-width) %>%
  dplyr::mutate(name = paste0(gene, "#", m_repeat))


relevant_snek_insertions <- sneks_in_exons %>%
  dplyr::mutate(name = paste0(gene, "#", m_repeat)) %>%
  dplyr::rename(start = repeat_start, end = repeat_end, strand = repeat_strand) %>%
  rbind(sneks_upstream_nongene)


relevant_snek_insertions_ranges <- plyranges::as_granges(relevant_snek_insertions)

# get sequence of upstream repeats to confirm
relevant_snek_insertions_seq <- Biostrings::getSeq(genome_seq, relevant_snek_insertions_ranges)
names(relevant_snek_insertions_seq) <- relevant_snek_insertions_ranges$name
Biostrings::writeXStringSet(relevant_snek_insertions_seq, "GeneInteraction/insertions.fa")

# peroform blast search against repeat annotation library
confirmation_blast_out <- read_tsv(system(paste0("blastn -dust yes -query GeneInteraction/insertions.fa -db ~/SnakeAnnotation/snake_annotation.lib -outfmt \"6 qseqid qstart qend sseqid sstart send pident qcovs bitscore length mismatch evalue qseqid slen qlen\""), intern = TRUE),
         col_names = c("gene", "qstart", "qend", "seqnames", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "slen", "qlen"))

# filter out repeats not identified as insertions
confirmed_repeats <- confirmation_blast_out %>%
  dplyr::mutate(gene_repeat = sub(".*#", "", gene), database_repeat = sub("#.*", "", seqnames)) %>%
  dplyr::filter(gene_repeat == database_repeat)

# create bed of insertions
confirmed_repeats_bed <- relevant_snek_insertions %>%
  dplyr::filter(name %in% confirmed_repeats$gene) %>%
  dplyr::mutate(X5 = ".") %>%
  dplyr::select(seqnames, start, end, name, X5, strand)

# write bed to file
readr::write_tsv(confirmed_repeats_bed, "GeneInteraction/speciesComparison/insertions.bed", col_names = F)
