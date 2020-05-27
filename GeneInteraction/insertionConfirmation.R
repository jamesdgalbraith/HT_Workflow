library(BSgenome)
library(GenomicRanges)
library(plyranges)
library(tidyverse)

# read bed of insertion locations
gene_list <- read_tsv("GeneInteraction/relevant_insertions_fixed.tsv") %>%
  mutate(appended_names = paste0(gene, "#", m_repeat), repeat_length = end - start + 1)

# read in genome sequence
genome_seq <- Biostrings::readDNAStringSet(filepath = "~/Genomes/Reptiles/Aipysurus_laevis/kmer_49.pilon_x2.sorted.fasta")
gc()
names(genome_seq) <- sub(" .*", "", names(genome_seq))
genome_idx <- tibble(seqnames = names(genome_seq), scaffold_length = as.double(width(genome_seq)))

# make ranges object of extended sequences
extended_repeat_ranges <- gene_list %>%
  as_tibble() %>%
  dplyr::mutate(start = start - 2000, end = end + 2000) %>%
  dplyr::inner_join(genome_idx) %>%
  dplyr::mutate(start = case_when(start < 1 ~ 1, start >= 1 ~ start),
         end = case_when(end > scaffold_length ~ scaffold_length, end <= scaffold_length ~ end)) %>%
  dplyr::select(-scaffold_length) %>%
  plyranges::as_granges()

# get sequence write fasta
extended_repeat_seq <- Biostrings::getSeq(genome_seq, extended_repeat_ranges)
names(extended_repeat_seq) <- extended_repeat_ranges$appended_names
Biostrings::writeXStringSet(extended_repeat_seq, "GeneInteraction/speciesComparison/temp.fa")

# search Emydocephalus
emydocephalus_blast_out <- read_tsv(system(paste0("blastn -dust yes -query GeneInteraction/speciesComparison/temp.fa -db ~/Genomes/Reptiles/Emydocephalus_ijimae/emyIji_1.0.fasta -outfmt \"6 qseqid qstart qend sseqid sstart send pident qcovs bitscore length mismatch evalue qseqid slen qlen\""), intern = TRUE),
                                    col_names = c("gene", "qstart", "qend", "seqnames", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "slen", "qlen")) 

# select best hits
emydocephalus_blast_out <- emydocephalus_blast_out %>%
  dplyr::filter(qcovs >= 40, pident >= 90) %>%
  dplyr::mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+"),
         start = case_when(sstart < send ~ sstart, send < sstart ~ send),
         end = case_when(sstart > send ~ sstart, send > sstart ~ send)) %>%
  dplyr::select(seqnames, start, end, strand, gene)

# convert into ranges object
emydocephalus_ranges <- plyranges::as_granges(emydocephalus_blast_out)

# reduce ranges object
emydocephalus_ranges_r <- GenomicRanges::reduce(emydocephalus_ranges, with.revmap = TRUE, min.gapwidth = 500)
revmap_val <- mcols(emydocephalus_ranges_r)$revmap
mcols(emydocephalus_ranges_r) <- DataFrame(gene = extractList(mcols(emydocephalus_ranges)$gene, revmap_val))

# select one hit per repeat
emydocephalus_cleaned <- emydocephalus_ranges_r %>%
  dplyr::mutate(gene = unlist(base::unique(gene))) %>%
  as_tibble() %>%
  group_by(gene) %>%
  dplyr::arrange(-width) %>%
  dplyr::slice(1) %>%
  ungroup()

# Hydrophis
hydrophis_blast_out <- read_tsv(system(paste0("blastn -dust yes -query GeneInteraction/speciesComparison/temp.fa -db ~/Genomes/Reptiles/Hydrophis_melanocephalus/hydMel_1.0.fasta -outfmt \"6 qseqid qstart qend sseqid sstart send pident qcovs bitscore length mismatch evalue qseqid\""), intern = TRUE),
                                col_names = c("gene", "qstart", "qend", "seqnames", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue")) 

# select best hits
hydrophis_blast_out <- hydrophis_blast_out %>%
  dplyr::filter(qcovs >= 40, pident >= 90) %>%
  dplyr::mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+"),
         start = case_when(sstart < send ~ sstart, send < sstart ~ send),
         end = case_when(sstart > send ~ sstart, send > sstart ~ send)) %>%
  dplyr::select(seqnames, start, end, strand, gene)

# convert into ranges object
hydrophis_ranges <- plyranges::as_granges(hydrophis_blast_out)

# reduce ranges object
hydrophis_ranges_r <- GenomicRanges::reduce(hydrophis_ranges, with.revmap = TRUE, min.gapwidth = 500)
revmap_val <- mcols(hydrophis_ranges_r)$revmap
mcols(hydrophis_ranges_r) <- DataFrame(gene = extractList(mcols(hydrophis_ranges)$gene, revmap_val))

# select one hit per repeat
hydrophis_cleaned <- hydrophis_ranges_r %>%
  dplyr::mutate(gene = unlist(base::unique(gene))) %>%
  as_tibble() %>%
  group_by(gene) %>%
  dplyr::arrange(-width) %>%
  dplyr::slice(1) %>%
  ungroup()

# Name genes and get sequences
aipysurus_cleaned <- extended_repeat_ranges %>%
  dplyr::mutate(gene = paste0("Aipysurus laevis ", appended_names))
aipysurus_cleaned_seq <- getSeq(genome_seq, aipysurus_cleaned)
names(aipysurus_cleaned_seq) <- aipysurus_cleaned$gene

# get emydocephalus sequence
genome_seq <- Biostrings::readDNAStringSet(filepath = "~/Genomes/Reptiles/Emydocephalus_ijimae/emyIji_1.0.fasta")
gc()
names(genome_seq) <- sub(" .*", "", names(genome_seq))

emydocephalus_cleaned <- emydocephalus_cleaned %>%
  dplyr::mutate(gene = paste0("Emydocephalus ijimae ", gene)) %>%
  plyranges::as_granges()
emydocephalus_cleaned_seq <- getSeq(genome_seq, emydocephalus_cleaned)
names(emydocephalus_cleaned_seq) <- emydocephalus_cleaned$gene

# get hydrophis sequence
genome_seq <- Biostrings::readDNAStringSet(filepath = "~/Genomes/Reptiles/Hydrophis_melanocephalus/hydMel_1.0.fasta")
gc()
names(genome_seq) <- sub(" .*", "", names(genome_seq))

hydrophis_cleaned <- hydrophis_cleaned %>%
  dplyr::mutate(gene = paste0("Hydrophis melanocephalus ", gene)) %>%
  plyranges::as_granges()
hydrophis_cleaned_seq <- getSeq(genome_seq, hydrophis_cleaned)
names(hydrophis_cleaned_seq) <- hydrophis_cleaned$gene

# compile sequences and align
for(i in 1:base::length(gene_list$name)){

  compiled_seq <- c(aipysurus_cleaned_seq[grepl(gene_list$appended_names[i], names(aipysurus_cleaned_seq))],
                    emydocephalus_cleaned_seq[grepl(gene_list$appended_names[i], names(emydocephalus_cleaned_seq))],
                    hydrophis_cleaned_seq[grepl(gene_list$appended_names[i], names(hydrophis_cleaned_seq))])
  
  if(base::length(compiled_seq) > 1){
    
    writeXStringSet(compiled_seq, "GeneInteraction/speciesComparison/temp.fa")
    system(paste0("mafft --adjustdirection --localpair --thread 12 GeneInteraction/speciesComparison/temp.fa > GeneInteraction/speciesComparison/", gene_list$appended_names[i], "_aligned.fa"))
    
  }
}

