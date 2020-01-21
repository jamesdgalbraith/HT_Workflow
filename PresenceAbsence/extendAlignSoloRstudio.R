# modified version of extendAlign.R for second stage (the extension)

suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

# set output folder
if(!dir.exists("curation/aligned")){dir.create("curation/aligned")}
out_folder <- "curation/"
genome_dir <- "~/Analysis/Genomes/"

# set variables
species_name <- "Notechis_scutatus"
genome_name <- "TS10Xv2-PRI.fasta"
print(species_name)

# set and read in genome and index
genome_path <- paste0(genome_dir, "/", species_name, "/", genome_name)

if(!file.exists(paste0(genome_path, ".fai"))){system(paste0("samtools faidx ", genome_path))}

genome_fai <- read_tsv(paste0(genome_path, ".fai"), col_names = c("seqnames", "scaffold_length", "x3", "x4", "x5")) %>%
  dplyr::select(1:2)

genome_seq <- Biostrings::readDNAStringSet(filepath = genome_path)
gc()
names(genome_seq) <- sub(" .*", "", names(genome_seq))

# set flanks, query and coverage as necessary
flank5 <- 1000
flank3 <- 1000
req_pident <- 90
req_length <- 500
query_repeat <- "../Aipysurus/Aipysurus_family117394_refined.fasta"

# set and read in repeats
repeat_seq <- Biostrings::readDNAStringSet(filepath = query_repeat)
names(repeat_seq) <- sub(" .*", "", names(repeat_seq))
repeat_name <- names(repeat_seq)

# For sequence from same species
blast_out <- read.table(text=system(paste0("blastn -query ", query_repeat, " -db ", genome_path, " -outfmt \"6 sseqid sstart send pident qcovs bitscore length\""), intern = TRUE), col.names = c("seqnames", "sstart", "send", "pident", "qcovs", "bitscore", "length"))

blast_out <- blast_out %>%
  as_tibble() %>%
  dplyr::mutate(seqnames = as.character(seqnames), sstart = as.double(sstart), send = as.double(send), length = as.double(length)) %>%
  dplyr::arrange(-bitscore) %>%
  mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+"),
         start = case_when(sstart < send ~ sstart, send < sstart ~ send),
         end = case_when(sstart > send ~ sstart, send > sstart ~ send))

# filter based on coverage and identity, extend and adjust
bed <- blast_out %>%
  inner_join(genome_fai) %>%
  dplyr::filter(length > req_length) %>%
  dplyr::filter(pident >= req_pident) %>%
  arrange(-pident) %>%
  mutate(start = case_when(strand == "+" ~ start - flank5, strand == "-" ~ start - flank3),
         end = case_when(strand == "-" ~ end + flank5, strand == "+" ~ end + flank3),
         start = case_when(start <= 1 ~ 1, start > 1 ~ start),
         end = case_when(end > scaffold_length ~ scaffold_length, end <= scaffold_length ~ end)) %>%
  dplyr::select(seqnames, start, end, strand)

# if more than 25 results select top 25
if(nrow(as_tibble(bed))>25){
  bed <- bed %>% 
    dplyr::slice(1:25)
}

bed <- plyranges::as_granges(bed)

bed <- reduce_ranges_directed(bed)

bed_tbl <- bed %>%
  as_tibble() %>%
  mutate(name = paste0(seqnames, ":", start, "-", end, "(", strand, ")"))


# get seqs
seqs <- Biostrings::getSeq(genome_seq, bed)

# name seqs
names(seqs) <- bed_tbl$name

# merge with original query
seqs <- c(repeat_seq, seqs)

# write seqs to files
Biostrings::writeXStringSet(x = seqs, filepath = paste0(out_folder, "/", species_name, "_", repeat_name, "_hits.fa"))

# perform multiple alignment
system(paste0("mafft --localpair --maxiterate 10 --thread 4 ", out_folder, "/", species_name, "_", repeat_name, "_hits.fa > ", out_folder, "/aligned/", species_name, "_", repeat_name, ".fa"))
