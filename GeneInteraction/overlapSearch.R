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
  inner_join(aipysurus_genome_annotation_seq_names)
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

mrna_overlap <- pair_overlaps(recip_ranges,mrna_gff_ranges) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames) %>%
  dplyr::rename(seqnames = granges.x.seqnames, repeat_start = granges.x.start, repeat_end = granges.x.end, repeat_strand = granges.x.strand,
                ann_start = granges.y.start, ann_end = granges.y.end, ann_strand = granges.y.strand) %>%
  mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand)) %>%
  dplyr::select(-granges.x.width, -granges.y.width) %>%
  base::unique()

exon_overlap <- pair_overlaps(recip_ranges,exon_gff_ranges) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames) %>%
  dplyr::rename(seqnames = granges.x.seqnames, repeat_start = granges.x.start, repeat_end = granges.x.end, repeat_strand = granges.x.strand,
                ann_start = granges.y.start, ann_end = granges.y.end, ann_strand = granges.y.strand) %>%
  mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand)) %>%
  dplyr::select(-granges.x.width, -granges.y.width) %>%
  base::unique()

mrna_overlap_repeat_bed <- mrna_overlap %>%
  dplyr::select(seqnames, repeat_start, repeat_end, repeat_strand, repeat_name) %>%
  dplyr::rename(start = repeat_start, end = repeat_end, strand = repeat_strand) %>%
  plyranges::as_granges() %>%
  plyranges::reduce_ranges_directed(repeat_name = base::unique(repeat_name)) %>%
  as_tibble() %>%
  mutate(repeat_name = unlist(repeat_name), X5 = ".") %>%
  select(seqnames, start, end, repeat_name, X5, strand)

write_tsv(mrna_overlap_repeat_bed, "GeneInteraction/mrna_overlap_repeat.bed", col_names = F)

mrna_overlap_annotation_bed <- mrna_gff %>%
  filter(gene %in% mrna_overlap$gene) %>%
  dplyr::mutate(X5 = ".") %>%
  dplyr::select(seqnames, start, end, gene, X5, strand)

write_tsv(mrna_overlap_annotation_bed, "GeneInteraction/mrna_overlap_annotation.bed", col_names = F)

exon_overlap_annotation_bed <- exon_gff %>%
  filter(gene %in% mrna_overlap$gene) %>%
  dplyr::mutate(X5 = ".") %>%
  dplyr::select(seqnames, start, end, gene, X5, strand)
  
write_tsv(exon_overlap_annotation_bed, "GeneInteraction/exon_overlap_annotation.bed", col_names = F)

mrna_overlap_repeat_bed %>%
  select(seqnames) %>%
  base::unique() %>%
  write_tsv("GeneInteraction/contigs_for_bams.txt", col_names = F)

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

utr_distance_gappy %>%
  filter(type == "five_prime_UTR") %>%
  filter((ann_strand == "+" & repeat_end < ann_start) | (ann_strand == "-" & ann_end < repeat_start)) %>%
  select(seqnames, repeat_start, repeat_end, repeat_strand, gene) %>%
  unique()

pair_overlaps(recip_ranges, aipysurus_genome_annotation) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames, -score, -phase, -X5_prime_partial, -X3_prime_partial, -Parent, -Name, -source) %>%
  rename(seqnames = granges.x.seqnames, repeat_start = granges.x.start, repeat_end = granges.x.end, ann_start = granges.y.start, ann_end = granges.y.end,
         repeat_strand = granges.x.strand, ann_strand = granges.y.strand) %>%
  mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand),
         ID = sub("\\.TU\\.", ".model.", ID), ID = sub("cds\\.", "", ID), ID = sub("\\.exon.*", "", ID), ID = sub(".utr.*", "", ID)) %>%
  inner_join(aipysurus_genome_annotation_seq_names) %>%
  select(-ID, -granges.y.width, -granges.x.width) %>%
  base::unique() %>%
  filter(grepl("CDS", type))

as_tibble(aipysurus_genome_annotation) %>%
  dplyr::select(-score, -phase, -X5_prime_partial, -X3_prime_partial, -Parent, -Name, -source) %>%
  dplyr::mutate(ID = sub("\\.TU\\.", ".model.", ID), ID = sub("cds\\.", "", ID), ID = sub("\\.exon.*", "", ID), ID = sub(".utr.*", "", ID)) %>%
  inner_join(aipysurus_genome_annotation_seq_names) %>%
  filter(grepl("CLIP4", gene)) %>%
  View()

# ACAA1 specific searches
ACAA1_aipysurus_ranges <- exon_overlap %>%
  dplyr::filter(gene == "ACAA1") %>%
  dplyr::mutate(start = repeat_start, end = repeat_end, strand = repeat_strand) %>%
  # dplyr::mutate(start = min(ann_start, repeat_start), end = max(ann_end, repeat_end), strand = ann_strand) %>%
  # dplyr::rename(start = ann_start, end = ann_end) %>%
  # dplyr::mutate(start = start - 2000, end = end + 2000) %>%
  plyranges::as_granges()

ACAA1_aipysurus_seq <- Biostrings::getSeq(genome_seq, ACAA1_aipysurus_ranges)
names(ACAA1_aipysurus_seq) <- paste0(seqnames(ACAA1_aipysurus_ranges), "-", ranges(ACAA1_aipysurus_ranges))
Biostrings::writeXStringSet(x = ACAA1_aipysurus_seq, filepath = "GeneInteraction/ACAA1_aipysurus_3utr.fa")


species_list <- read_tsv("snake_species.tsv", col_names = c("s_name", "g_name"))

# loop to find other hits
for (j in 1:nrow(species_list)){
  # print name of species
  print(species_list$s_name[j])
  if(j == 1){
    compiled_ACAA1_seq <- ACAA1_aipysurus_seq
    next()
  }
  # run blast, modified from normal to cope with one line results
  other_blast_out <- system(paste0("blastn -dust yes -num_threads 12 -query GeneInteraction/ACAA1_aipysurus_3utr.fa -db ~/Genomes/Reptiles/", species_list$s_name[j], "/", species_list$g_name[j], " -outfmt \"6 qseqid sseqid qstart qend sstart send pident qcovs bitscore length mismatch evalue\""), intern = TRUE) %>%
    tibble() %>%
    tidyr::separate(1, sep = "\t",
                    into = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue"))
  
  # read and index genome
  other_genome_seq <- Biostrings::readDNAStringSet(filepath = paste0(genome_dir, "/", species_list$s_name[j], "/", species_list$g_name[j]))
  gc()
  names(other_genome_seq) <- sub(" .*", "", names(other_genome_seq))
  other_genome_idx <- tibble(seqnames = names(other_genome_seq), width = as.double(width(other_genome_seq)))
  
  # convert best blast hit to ranges
  ACAA1_other_ranges <- other_blast_out %>%
    dplyr::slice(1:2) %>%
    dplyr::mutate(sstart = as.integer(sstart), send = as.integer(send), qstart = as.integer(qstart), qend = as.integer(qend),
                  strand = case_when(sstart < send ~ "+", sstart > send ~ "-"),
                  start = case_when(strand == "+" ~ sstart, strand == "-" ~ send),
                  end = case_when(strand == "+" ~ send, strand == "-" ~ sstart)) %>%
    dplyr::rename(seqnames = sseqid) %>%
    # dplyr::mutate(start = start - 2000, end = end + 2000) %>%
    # inner_join(other_genome_idx) %>%
    # dplyr::mutate(start = case_when(start >= 1 ~ start, start < 1 ~ 1),
    #               end = case_when(end <= width ~ end, end > width ~ width)) %>%
    dplyr::select(seqnames, start, end, strand) %>%
    plyranges::as_granges() %>%
    IRanges::reduce(min.gapwidth=2000)
  
  ACAA1_other_seq <- Biostrings::getSeq(other_genome_seq, ACAA1_other_ranges)
  names(ACAA1_other_seq) <- paste0(species_list$s_name[j], "_", seqnames(ACAA1_other_ranges), "-", ranges(ACAA1_other_ranges))
  compiled_ACAA1_seq <- c(compiled_ACAA1_seq, ACAA1_other_seq)
 
}

Biostrings::writeXStringSet(x = compiled_ACAA1_seq, filepath = "GeneInteraction/ACAA1_elpaid_extended_3utr.fa")

system("mafft --thread 12 GeneInteraction/ACAA1_elpaid_extended_3utr.fa > GeneInteraction/ACAA1_elpaid_extended_3utr_aligned.fa")
names(compiled_ACAA1_seq)

Biostrings::writeXStringSet(x = compiled_ACAA1_seq, filepath = "GeneInteraction/ACAA1_elpaid_extended_search_3utr.fa")

system("mafft --thread 12 GeneInteraction/ACAA1_elpaid_extended_search_3utr.fa > GeneInteraction/ACAA1_elpaid_extended_search_3utr_aligned.fa")
names(compiled_ACAA1_seq)
