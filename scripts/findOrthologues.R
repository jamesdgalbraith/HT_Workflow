suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

# set variables
query_path <- "~/Analysis/Snake/HT_Workflow/RTE-Snek/RTE-Snek.fasta"
primary_genome <- "~/Analysis/Genomes/Aipysurus_laevis/assembly_20171114.fasta"

# query_path <- args$query
# genome <- args$genome
# repeat_database <- args$repeat_database

# create names of query and subject based on query path
query_name <- sub(".*\\/", "", query_path)
query_name <- sub("\\.fasta", "", query_name)
primary_genome_name <- sub(".*\\/", "", primary_genome)
primary_genome_name <- sub("\\.fasta", "", primary_genome_name)

# ensure blast db of genome exists
if(!file.exists(paste0(primary_genome, ".nhr")) | !file.exists(paste0(primary_genome, ".nin")) | !file.exists(paste0(primary_genome, ".nsq"))){
  system(paste0("makeblastdb -dbtype nucl -in ", primary_genome, " -out ", primary_genome))
}

# ensure index of genome exists
if(!file.exists(paste0(primary_genome, ".fai"))){
  system(paste0("samtools faidx ", primary_genome))
}

# read in index of genome
genome_fai <- readr::read_tsv(paste0(primary_genome, ".fai"), col_names = c("sseqid", "scaffold_length", "X1", "X2", "X3")) %>% select(1:2)

# ensure working directory exists
if(!dir.exists("./working")){dir.create("./working")}

# run blast / read in blast table
system(paste0("blastn -evalue 1.0e-100 -num_threads 4 -outfmt 6 -dust yes -query ", query_path, " -db ", primary_genome, " -out working/", query_name, "_in_", primary_genome_name, ".out"))

blast_1 <- read_tsv(paste0("working/", query_name, "_in_", primary_genome_name, ".out"), col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")) %>%
  inner_join(genome_fai)

# manipulate bl;ast data into useable format
blast_1 <- blast_1 %>%
  dplyr::mutate(strand = base::ifelse(sstart < send, "+", "-")) %>%
  dplyr::mutate(start = base::ifelse(sstart < send, sstart, send)) %>%
  dplyr::mutate(end = base::ifelse(sstart > send, send, sstart)) %>%
  dplyr::rename(seqnames = sseqid)

# filter out small hits
blast_1 <- blast_1 %>%
  filter(length >= 1000)

# extend queries
blast_ext <- blast_1 %>%
  mutate(start_ext = base::ifelse( start <= 10000, 1, sstart - 10000),
         end_ext = base::ifelse(end + 10000 >= scaffold_length, scaffold_length, end + 10000)) %>%
  dplyr::mutate(names = paste0(primary_genome_name, "#", seqnames, ":", start_ext, "-", end_ext, "(", strand, ")"))

# convert extended to ranges
blast_ext_ranges <- blast_ext %>%
  dplyr::select(seqnames, start_ext, end_ext, strand, names) %>%
  dplyr::rename(start = start_ext, end = end_ext) %>%
  plyranges::as_granges()

# read in genome, rename seq, collect garbage
genome_seq <- Biostrings::readDNAStringSet(filepath = primary_genome)
names(genome_seq) <- sub(" .*", "", names(genome_seq))
gc(verbose = F)

# get sequences
blast_1_seqs <- Biostrings::getSeq(genome_seq, blast_ext_ranges)

# names sequences
names(blast_1_seqs) <- blast_ext_ranges$names

blast_ext_ranges %>%
  as_tibble %>%
  mutate(qseqid = names, genome_name = primary_genome_name, names = paste0(genome_name, "#", names))

blast_ext_ranges$seq <- blast_1_seqs

# write sequences to file
Biostrings::writeXStringSet(x = blast_1_seqs, filepath = paste0("working/", query_name, "_in_", primary_genome_name, ".fasta"))

#### Search other genomes
comparison_genomes <- readr::read_tsv("~/Analysis/Snake/HT_Workflow/comparison_genomes.txt", col_names = "genome_path")

# subject genome
secondary_genome <- comparison_genomes$genome_path[2]
subject_name <- sub(".*\\/", "", secondary_genome)
subject_name <- sub("\\.fasta", "", subject_name)

# ensure blast db of genome exists
if(!file.exists(paste0(secondary_genome, ".nhr")) | !file.exists(paste0(secondary_genome, ".nin")) | !file.exists(paste0(secondary_genome, ".nsq"))){
  system(paste0("makeblastdb -dbtype nucl -in ", secondary_genome, " -out ", secondary_genome))
}

# ensure index of genome exists
if(!file.exists(paste0(secondary_genome, ".fai"))){
  system(paste0("samtools faidx ", secondary_genome))
}

# search - blast sequences against subject genome
system(paste0("blastn -query working/", query_name, "_in_", primary_genome_name, ".fasta -db ", secondary_genome, " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs\" -out working/", query_name, "_in_", subject_name, ".out"))

# read in reciprocal search                      
blast_2 <- read_tsv(file = paste0("working/", query_name, "_in_", subject_name, ".out"), col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs")) %>%
  dplyr::mutate(strand = base::ifelse(sstart < send, "+", "-")) %>%
  dplyr::mutate(start = base::ifelse(strand == "+", sstart, send)) %>%
  dplyr::mutate(end = base::ifelse(strand == "+", send, sstart)) %>%
  dplyr::select(-sstart, -send)

# filter out hits with low coverage, identity and hit length
blast_2 <- blast_2 %>%
  dplyr::filter(qcovs > 30, length > 1000, pident >= 90) %>%
  dplyr::arrange(qseqid, sseqid, start)

# select smallest and largest hit on genome
blast_2_ranges <- blast_2 %>%
  dplyr::group_by(qseqid, sseqid, strand) %>%
  dplyr::summarise(start = min(start), end = max(end)) %>%
  dplyr::mutate(genome_name = subject_name, names = paste0(sseqid, ":", start, "-", end, "(", strand, ")")) %>%
  dplyr::rename(seqnames = sseqid) %>%
  plyranges::as_granges()

# read in second genome
genome_seq <- Biostrings::readDNAStringSet(filepath = secondary_genome)
names(genome_seq) <- sub(" .*", "", names(genome_seq))
gc(verbose = F)

# get subject ranges
blast_2_seqs <- Biostrings::getSeq(genome_seq, blast_2_ranges)

# name sequences
names(blast_2_seqs) <- blast_2_ranges$names

# add seq to ranges object
blast_2_ranges$seq <- blast_2_seqs

# finalise object for output
blast_2_final <- blast_2_ranges %>%
  as_tibble()


blast_ext_ranges %>%
  as_tibble %>%
  mutate(qseqid = names, genome_name = primary_genome_name, names = paste0(genome_name, "#", names), seqs = blast)

