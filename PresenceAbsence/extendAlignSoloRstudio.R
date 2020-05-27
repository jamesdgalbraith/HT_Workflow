# modified version of extendAlign.R for second stage (the extension)

suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

# query repeat family name
query_repeat_family_name <- "RTE-Snek_2"

# set output folder
out_folder <- paste0("../", query_repeat_family_name, "/curation_step2/aligned")
if(!dir.exists(out_folder)){dir.create(out_folder, recursive = T)}

# set variables
# set and read in genome and index
genome_path <- "~/Genomes/Reptiles/Hydrophis_melanocephalus/hydMel_1.0.fasta"
genome_seq <- Biostrings::readDNAStringSet(filepath = genome_path)
gc()
names(genome_seq) <- sub(" .*", "", names(genome_seq))
genome_fai <- tibble(seqnames = names(genome_seq), scaffold_length = as.double(width(genome_seq)))


# set flanks, query and coverage as necessary
flank5 <- 0
flank3 <- 0
req_pident <- 96
req_length <- 50

query_repeat <- "~/HT_Workflow/RTE-Snek_2/curation/curated_1/Hydrophis_melanocephalus_RTE-Snek_2.fasta"

# set and read in repeats
repeat_seq <- Biostrings::readDNAStringSet(filepath = query_repeat)
names(repeat_seq) <- sub(" .*", "", names(repeat_seq))
repeat_name <- names(repeat_seq)

# For sequence from same species
blast_out <- read.table(text=system(paste0("blastn -query ", query_repeat, " -db ", genome_path, " -outfmt \"6 sseqid sstart send pident qcovs bitscore length qstart qend\""), intern = TRUE), col.names = c("seqnames", "sstart", "send", "pident", "qcovs", "bitscore", "length", "qstart", "qend"))

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
  dplyr::filter(pident >= req_pident) %>%
  dplyr::filter(length >= req_length) %>%
  arrange(-bitscore) %>%
  mutate(start = case_when(strand == "+" ~ start - flank5, strand == "-" ~ start - flank3),
         end = case_when(strand == "-" ~ end + flank5, strand == "+" ~ end + flank3)) %>%
  mutate(start = case_when(start <= 1 ~ 1, start > 1 ~ start),
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
Biostrings::writeXStringSet(x = seqs, filepath = paste0(out_folder, "/", repeat_name, "_hits.fa"))

# perform multiple alignment
system(paste0("mafft --localpair --maxiterate 10 --thread 12 ", out_folder, "/", repeat_name, "_hits.fa > ", out_folder, "/", repeat_name, "_aligned.fa"))

