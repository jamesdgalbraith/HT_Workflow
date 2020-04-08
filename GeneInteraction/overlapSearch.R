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

aipysurus_sneks_annotation <- aipysurus_repeat_annotation_corrected %>%
  filter(grepl("Snek", full_repeat), end - start > 50, perc_div <= 25)

aipysurus_sneks_annotation_ranges <- plyranges::as_granges(aipysurus_sneks_annotation)

aipysurus_lines_annotation_ranges <- plyranges::as_granges(aipysurus_lines_annotation)

flattened <- aipysurus_lines_annotation_ranges %>% GenomicRanges::reduce(min.gapwidth = 500)

tibble(name = names(table(aipysurus_lines_annotation$full_repeat)), count = table(aipysurus_lines_annotation$full_repeat)) %>%
  filter(count > 10) %>%
  arrange(-count)

ggplot(aipysurus_lines_annotation, aes(x = (end - start + 1), y = perc_div)) + geom_point()
ggplot(aipysurus_lines_annotation, aes(x = perc_div)) + geom_histogram(bins = 30)
ggplot(aipysurus_lines_annotation, aes(x = (end - start + 1))) + geom_histogram()


# # Search for lines in genome
# blast_out <- read_tsv(system(paste0("blastn -dust yes -query compiled_lines.fasta -db ", genome_path, " -outfmt \"6 sseqid sstart send pident qcovs bitscore length mismatch evalue qseqid\""), intern = TRUE), col_names = c("seqnames", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "qseqid")) %>%
#   filter(length > 50, pident >= 97) %>%
#   mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+"),
#          start = case_when(sstart < send ~ sstart, send < sstart ~ send),
#          end = case_when(sstart > send ~ sstart, send > sstart ~ send))
# 
# # Create bed of blast output
# bed <- blast_out %>%
#   inner_join(genome_idx) %>%
#   dplyr::select(seqnames, start, end, strand) %>%
#   as_granges() %>%
#   GenomicRanges::reduce(min.gapwidth = 500)
# bed$names <- paste0(seqnames(bed), ":", ranges(bed), "(", strand(bed), ")")
# 
# # Get seqs of blast hits and write to file
# seqs <- Biostrings::getSeq(genome_seq, bed)
# names(seqs) <- bed$names
# Biostrings::writeXStringSet(x = seqs, filepath = "GeneInteraction/initialHits.fa")
# 
# # align hits to original search and select best hit for each based on bitscore
# blastn_recip <- read_tsv(system(paste0("blastn -dust yes -query GeneInteraction/initialHits.fa -subject compiled_lines.fasta  -outfmt \"6 qseqid sseqid qstart qend sstart send pident qcovs bitscore length mismatch evalue\""), intern = TRUE), col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue")) %>%
#   group_by(qseqid) %>%
#   dplyr::slice(1) %>%
#   select(qseqid, sseqid) %>%
#   dplyr::rename(names = qseqid, repeat_name = sseqid) %>%
#   ungroup()
# recip_tibble <- bed %>%
#   as_tibble() %>%
#   mutate(seqnames = as.character(seqnames)) %>%
#   dplyr::rename(hit_width = width) %>%
#   inner_join(blastn_recip)
# recip_ranges <- recip_tibble %>%
#   select(-hit_width, -names) %>%
#   as_granges()

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

nearest_line_gene_rev <- pair_precede(gene_gff_ranges, aipysurus_lines_annotation_ranges) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames) %>%
  rename(seqnames = granges.x.seqnames, ann_start = granges.x.start, ann_end = granges.x.end, repeat_start = granges.y.start, repeat_end = granges.y.end,
         ann_strand = granges.x.strand, repeat_strand = granges.y.strand) %>%
  mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand),
         dist = repeat_start - ann_end + 1) %>%
  filter(ann_strand == "-") %>%
  base::unique()

nearest_line_gene_fwd <- pair_follow(gene_gff_ranges, aipysurus_lines_annotation_ranges) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames) %>%
  rename(seqnames = granges.x.seqnames, ann_start = granges.x.start, ann_end = granges.x.end, repeat_start = granges.y.start, repeat_end = granges.y.end,
         ann_strand = granges.x.strand, repeat_strand = granges.y.strand) %>%
  mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand),
         dist = ann_start - repeat_end + 1) %>%
  filter(ann_strand == "+") %>%
  base::unique()

nearest_line_gene <- base::rbind(nearest_line_gene_fwd, nearest_line_gene_rev) %>%
  dplyr::arrange(seqnames, ann_start)

nearest_line_gene %>%
  filter(grepl("Snek", full_repeat))

nearest_gene_to_lines <- pair_overlaps(gene_gff_ranges, aipysurus_lines_annotation_ranges, maxgap = 5000) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames) %>%
  rename(seqnames = granges.x.seqnames, ann_start = granges.y.start, ann_end = granges.y.end, repeat_start = granges.x.start, repeat_end = granges.x.end,
         ann_strand = granges.y.strand, repeat_strand = granges.x.strand) %>%
  mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand),
         dist = case_when(ann_strand == "+" ~ , ann_strand == "-")) %>%
  base::unique()

lines_in_genes <- pair_overlaps(gene_gff_ranges, aipysurus_lines_annotation_ranges) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames) %>%
  rename(seqnames = granges.x.seqnames, ann_start = granges.y.start, ann_end = granges.y.end, repeat_start = granges.x.start, repeat_end = granges.x.end,
         ann_strand = granges.y.strand, repeat_strand = granges.x.strand) %>%
  mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand)) %>%
  

pair_overlaps(aipysurus_lines_annotation_ranges, exon_gff_ranges) %>%
  as_tibble() %>%
  dplyr::select(-granges.y.seqnames) %>%
  rename(seqnames = granges.x.seqnames, ann_start = granges.y.start, ann_end = granges.y.end, repeat_start = granges.x.start, repeat_end = granges.x.end,
         ann_strand = granges.y.strand, repeat_strand = granges.x.strand) %>%
  mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand))


# mrna_overlap <- pair_overlaps(recip_ranges,mrna_gff_ranges) %>%
#   as_tibble() %>%
#   dplyr::select(-granges.y.seqnames) %>%
#   dplyr::rename(seqnames = granges.x.seqnames, repeat_start = granges.x.start, repeat_end = granges.x.end, repeat_strand = granges.x.strand,
#                 ann_start = granges.y.start, ann_end = granges.y.end, ann_strand = granges.y.strand) %>%
#   mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand)) %>%
#   dplyr::select(-granges.x.width, -granges.y.width) %>%
#   base::unique()
# 
# exon_overlap <- pair_overlaps(recip_ranges,exon_gff_ranges) %>%
#   as_tibble() %>%
#   dplyr::select(-granges.y.seqnames) %>%
#   dplyr::rename(seqnames = granges.x.seqnames, repeat_start = granges.x.start, repeat_end = granges.x.end, repeat_strand = granges.x.strand,
#                 ann_start = granges.y.start, ann_end = granges.y.end, ann_strand = granges.y.strand) %>%
#   mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand)) %>%
#   dplyr::select(-granges.x.width, -granges.y.width) %>%
#   base::unique()
# 
# mrna_overlap_repeat_bed <- mrna_overlap %>%
#   dplyr::select(seqnames, repeat_start, repeat_end, repeat_strand, repeat_name) %>%
#   dplyr::rename(start = repeat_start, end = repeat_end, strand = repeat_strand) %>%
#   plyranges::as_granges() %>%
#   plyranges::reduce_ranges_directed(repeat_name = base::unique(repeat_name)) %>%
#   as_tibble() %>%
#   mutate(repeat_name = unlist(repeat_name), X5 = ".") %>%
#   select(seqnames, start, end, repeat_name, X5, strand)
# 
# write_tsv(mrna_overlap_repeat_bed, "GeneInteraction/mrna_overlap_repeat.bed", col_names = F)
# 
# mrna_overlap_annotation_bed <- mrna_gff %>%
#   filter(gene %in% mrna_overlap$gene) %>%
#   dplyr::mutate(X5 = ".") %>%
#   dplyr::select(seqnames, start, end, gene, X5, strand)
# 
# write_tsv(mrna_overlap_annotation_bed, "GeneInteraction/mrna_overlap_annotation.bed", col_names = F)
# 
# exon_overlap_annotation_bed <- exon_gff %>%
#   filter(gene %in% mrna_overlap$gene) %>%
#   dplyr::mutate(X5 = ".") %>%
#   dplyr::select(seqnames, start, end, gene, X5, strand)
#   
# write_tsv(exon_overlap_annotation_bed, "GeneInteraction/exon_overlap_annotation.bed", col_names = F)
# 
# mrna_overlap_repeat_bed %>%
#   select(seqnames) %>%
#   base::unique() %>%
#   write_tsv("GeneInteraction/contigs_for_bams.txt", col_names = F)
# 
# utr_distance_gappy <- pair_overlaps(recip_ranges, aipysurus_genome_annotation, maxgap = 5000) %>%
#   as_tibble() %>%
#   dplyr::select(-granges.y.seqnames, -score, -phase, -X5_prime_partial, -X3_prime_partial, -Parent, -Name, -source) %>%
#   rename(seqnames = granges.x.seqnames, repeat_start = granges.x.start, repeat_end = granges.x.end, ann_start = granges.y.start, ann_end = granges.y.end,
#          repeat_strand = granges.x.strand, ann_strand = granges.y.strand) %>%
#   mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand),
#          ID = sub("\\.TU\\.", ".model.", ID), ID = sub("cds\\.", "", ID), ID = sub("\\.exon.*", "", ID), ID = sub(".utr.*", "", ID)) %>%
#   inner_join(aipysurus_genome_annotation_seq_names) %>%
#   select(-ID, -granges.y.width, -granges.x.width) %>%
#   base::unique() %>%
#   filter(grepl("UTR", type))
# 
# utr_distance_gappy %>%
#   filter(type == "five_prime_UTR") %>%
#   filter((ann_strand == "+" & repeat_end < ann_start) | (ann_strand == "-" & ann_end < repeat_start)) %>%
#   select(seqnames, repeat_start, repeat_end, repeat_strand, gene) %>%
#   unique()
# 
# pair_overlaps(recip_ranges, aipysurus_genome_annotation) %>%
#   as_tibble() %>%
#   dplyr::select(-granges.y.seqnames, -score, -phase, -X5_prime_partial, -X3_prime_partial, -Parent, -Name, -source) %>%
#   rename(seqnames = granges.x.seqnames, repeat_start = granges.x.start, repeat_end = granges.x.end, ann_start = granges.y.start, ann_end = granges.y.end,
#          repeat_strand = granges.x.strand, ann_strand = granges.y.strand) %>%
#   mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand),
#          ID = sub("\\.TU\\.", ".model.", ID), ID = sub("cds\\.", "", ID), ID = sub("\\.exon.*", "", ID), ID = sub(".utr.*", "", ID)) %>%
#   inner_join(aipysurus_genome_annotation_seq_names) %>%
#   select(-ID, -granges.y.width, -granges.x.width) %>%
#   base::unique() %>%
#   filter(grepl("CDS", type))
# 
# ##### other repeats #####
# aipysurus_lines_annotation_sneks <- aipysurus_lines_annotation %>%
#   filter(grepl("Snek", full_repeat)) %>%
#   dplyr::select(seqnames, start, end, strand, full_repeat, perc_div)
# 
# aipysurus_lines_annotation_snek_ranges <- aipysurus_lines_annotation_sneks %>%
#   plyranges::as_granges()
# 
# aipysurus_lines_annotation_snekless <- aipysurus_lines_annotation %>%
#   filter(!grepl("Snek", full_repeat)) %>%
#   dplyr::select(seqnames, start, end, strand, full_repeat, perc_div)
# 
# aipysurus_lines_annotation_snekless_ranges <- aipysurus_lines_annotation_snekless %>%
#   plyranges::as_granges()
# 
# snekless_overlaps <- pair_overlaps(aipysurus_lines_annotation_snekless_ranges, exon_gff_ranges, maxgap = 5000) %>%
#   as_tibble() %>%
#   dplyr::select(-granges.y.seqnames) %>%
#   rename(seqnames = granges.x.seqnames, repeat_start = granges.x.start, repeat_end = granges.x.end, ann_start = granges.y.start, ann_end = granges.y.end,
#          repeat_strand = granges.x.strand, ann_strand = granges.y.strand, repeat_width = granges.x.width, ann_width = granges.y.width) %>%
#   mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand))  %>%
#   base::unique()
# 
# snek_overlaps <- pair_overlaps(aipysurus_lines_annotation_snek_ranges, exon_gff_ranges) %>%
#   as_tibble() %>%
#   dplyr::select(-granges.y.seqnames) %>%
#   rename(seqnames = granges.x.seqnames, repeat_start = granges.x.start, repeat_end = granges.x.end, ann_start = granges.y.start, ann_end = granges.y.end,
#          repeat_strand = granges.x.strand, ann_strand = granges.y.strand, repeat_width = granges.x.width, ann_width = granges.y.width) %>%
#   mutate(seqnames = as.character(seqnames), ann_strand = as.character(ann_strand), repeat_strand = as.character(repeat_strand))  %>%
#   base::unique()
# 
# snekless_overlaps_reduced <- snekless_overlaps %>%
#   dplyr::rename(start = repeat_start, end = repeat_end, strand = repeat_strand) %>%
#   dplyr::select(seqnames, start, end, strand) %>%
#   as_granges() %>%
#   GenomicRanges::reduce(min.gapwidth = 100, ignore.strand=F)
# 
# snek_overlaps_reduced <- snek_overlaps %>%
#   dplyr::rename(start = repeat_start, end = repeat_end, strand = repeat_strand) %>%
#   dplyr::select(seqnames, start, end, strand) %>%
#   as_granges() %>%
#   GenomicRanges::reduce(min.gapwidth = 100, ignore.strand=F) %>%
#   as_tibble()
# 
# snekless_overlaps_reduced_seq <- Biostrings::getSeq(genome_seq, snekless_overlaps_reduced)
# 
# names(snekless_overlaps_reduced_seq) <- paste0(seqnames(snekless_overlaps_reduced), ":", ranges(snekless_overlaps_reduced), "(", strand(snekless_overlaps_reduced), ")")
# 
# Biostrings::writeXStringSet(snekless_overlaps_reduced_seq, "GeneInteraction/snekless_repeats.fa")
# 
# 
# 
# ##### specific gene interactions, (ignore as necessary) #####
# # as_tibble(aipysurus_genome_annotation) %>%
# #   dplyr::select(-score, -phase, -X5_prime_partial, -X3_prime_partial, -Parent, -Name, -source) %>%
# #   dplyr::mutate(ID = sub("\\.TU\\.", ".model.", ID), ID = sub("cds\\.", "", ID), ID = sub("\\.exon.*", "", ID), ID = sub(".utr.*", "", ID)) %>%
# #   inner_join(aipysurus_genome_annotation_seq_names) %>%
# #   filter(grepl("CLIP4", gene)) %>%
# #   View()
# # 
# # # ACAA1 specific searches
# # ACAA1_aipysurus_ranges <- exon_overlap %>%
# #   dplyr::filter(gene == "ACAA1") %>%
# #   dplyr::mutate(start = repeat_start, end = repeat_end, strand = repeat_strand) %>%
# #   # dplyr::mutate(start = min(ann_start, repeat_start), end = max(ann_end, repeat_end), strand = ann_strand) %>%
# #   # dplyr::rename(start = ann_start, end = ann_end) %>%
# #   # dplyr::mutate(start = start - 2000, end = end + 2000) %>%
# #   plyranges::as_granges()
# # 
# # ACAA1_aipysurus_seq <- Biostrings::getSeq(genome_seq, ACAA1_aipysurus_ranges)
# # names(ACAA1_aipysurus_seq) <- paste0(seqnames(ACAA1_aipysurus_ranges), "-", ranges(ACAA1_aipysurus_ranges))
# # Biostrings::writeXStringSet(x = ACAA1_aipysurus_seq, filepath = "GeneInteraction/ACAA1_aipysurus_3utr.fa")
# # 
# # 
# # species_list <- read_tsv("snake_species.tsv", col_names = c("s_name", "g_name"))
# # 
# # # loop to find other hits
# # for (j in 1:nrow(species_list)){
# #   # print name of species
# #   print(species_list$s_name[j])
# #   if(j == 1){
# #     compiled_ACAA1_seq <- ACAA1_aipysurus_seq
# #     next()
# #   }
# #   # run blast, modified from normal to cope with one line results
# #   other_blast_out <- system(paste0("blastn -dust yes -num_threads 12 -query GeneInteraction/ACAA1_aipysurus_3utr.fa -db ~/Genomes/Reptiles/", species_list$s_name[j], "/", species_list$g_name[j], " -outfmt \"6 qseqid sseqid qstart qend sstart send pident qcovs bitscore length mismatch evalue\""), intern = TRUE) %>%
# #     tibble() %>%
# #     tidyr::separate(1, sep = "\t",
# #                     into = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue"))
# #   
# #   # read and index genome
# #   other_genome_seq <- Biostrings::readDNAStringSet(filepath = paste0(genome_dir, "/", species_list$s_name[j], "/", species_list$g_name[j]))
# #   gc()
# #   names(other_genome_seq) <- sub(" .*", "", names(other_genome_seq))
# #   other_genome_idx <- tibble(seqnames = names(other_genome_seq), width = as.double(width(other_genome_seq)))
# #   
# #   # convert best blast hit to ranges
# #   ACAA1_other_ranges <- other_blast_out %>%
# #     dplyr::slice(1:2) %>%
# #     dplyr::mutate(sstart = as.integer(sstart), send = as.integer(send), qstart = as.integer(qstart), qend = as.integer(qend),
# #                   strand = case_when(sstart < send ~ "+", sstart > send ~ "-"),
# #                   start = case_when(strand == "+" ~ sstart, strand == "-" ~ send),
# #                   end = case_when(strand == "+" ~ send, strand == "-" ~ sstart)) %>%
# #     dplyr::rename(seqnames = sseqid) %>%
# #     # dplyr::mutate(start = start - 2000, end = end + 2000) %>%
# #     # inner_join(other_genome_idx) %>%
# #     # dplyr::mutate(start = case_when(start >= 1 ~ start, start < 1 ~ 1),
# #     #               end = case_when(end <= width ~ end, end > width ~ width)) %>%
# #     dplyr::select(seqnames, start, end, strand) %>%
# #     plyranges::as_granges() %>%
# #     IRanges::reduce(min.gapwidth=2000)
# #   
# #   ACAA1_other_seq <- Biostrings::getSeq(other_genome_seq, ACAA1_other_ranges)
# #   names(ACAA1_other_seq) <- paste0(species_list$s_name[j], "_", seqnames(ACAA1_other_ranges), "-", ranges(ACAA1_other_ranges))
# #   compiled_ACAA1_seq <- c(compiled_ACAA1_seq, ACAA1_other_seq)
# #  
# # }
# # 
# # Biostrings::writeXStringSet(x = compiled_ACAA1_seq, filepath = "GeneInteraction/ACAA1_elpaid_extended_3utr.fa")
# # 
# # system("mafft --thread 12 GeneInteraction/ACAA1_elpaid_extended_3utr.fa > GeneInteraction/ACAA1_elpaid_extended_3utr_aligned.fa")
# # names(compiled_ACAA1_seq)
# # 
# # Biostrings::writeXStringSet(x = compiled_ACAA1_seq, filepath = "GeneInteraction/ACAA1_elpaid_extended_search_3utr.fa")
# # 
# # system("mafft --thread 12 GeneInteraction/ACAA1_elpaid_extended_search_3utr.fa > GeneInteraction/ACAA1_elpaid_extended_search_3utr_aligned.fa")
# # names(compiled_ACAA1_seq)
