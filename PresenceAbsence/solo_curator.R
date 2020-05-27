# read in libraries, hiding messages
suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

# set species and genome names
clade_name <- "Reptiles"
species_name <- "Aipysurus_laevis"
genome_name <- "ASM378908v1.fasta"

# set/create output directory
out_dir <- paste0("PresenceAbsence/curation/RTE-Snek_2/")
if(!dir.exists(out_dir))dir.create(out_dir)

# set genome and carp paths
genome_path <- paste0("~/Genomes/", clade_name, "/", species_name, "/", genome_name)

# read in genome
genome_seq <- Biostrings::readDNAStringSet(filepath = genome_path)
gc()
names(genome_seq) <- sub(pattern = " .*", replacement = "", x = names(genome_seq))

# create tibble of genome sequence names and lengths
genome_fai <- tibble(seqnames = names(genome_seq), scaffold_length = as.double(width(genome_seq)))

# set query
repeat_name <- "Penaeus_vannamei_RTE-Snek_2b.fasta"
query_repeat <- paste0("PresenceAbsence/curation/RTE-Snek_2/", repeat_name)
query_seq <- readDNAStringSet(query_repeat)
# set flanks either side
flank5 <- 100
flank3 <- 100
  
# blast search for repeat
blast_out <- read_tsv(system(paste0("blastn -num_threads 12 -query ", query_repeat, " -db ", genome_path, " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs\""), intern = T), col_names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs")) %>%
    dplyr::filter(length >= 0.25 * qlen) %>%
    dplyr::mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+"), # determine strand of hit
                  start = case_when(strand == "+" ~ sstart, strand == "-" ~ send), # determine stranded start
                  end = case_when(strand == "+" ~ send, strand == "-" ~ sstart))

  # subest base on query repeat
  query_blast_subset <- blast_out %>%
    dplyr::arrange(-bitscore)
  
  blast_out
  
  # reduce number of hits if > 20
  if(base::nrow(query_blast_subset) > 20){
    query_blast_subset <- query_blast_subset %>%
      dplyr::slice(1:20)
  }
  
  # adjust and create granges object
  query_blast_subset_ranges <- query_blast_subset %>%
    dplyr::select(seqnames, start, end, strand, slen) %>%
    dplyr::mutate(start = case_when(strand == "+" ~ start - flank5, strand == "-" ~ start - flank3),
                  end = case_when(strand == "+" ~ end + flank3, strand == "-" ~ end + flank5),
                  start = case_when(start > 0 ~ start, start <= 0 ~ 1),
                  end = case_when(end <= slen ~ end, end > slen ~ slen)) %>%
    plyranges::as_granges()
  
  query_blast_subset_ranges <- reduce_ranges_directed(query_blast_subset_ranges)
  
  # get sequences from genome
  query_blast_subset_seq <- Biostrings::getSeq(genome_seq, query_blast_subset_ranges)
  
  # name sequences using scaffold_name:cordinates(strand) format
  names(query_blast_subset_seq) <- paste0(seqnames(query_blast_subset_ranges), ":", ranges(query_blast_subset_ranges), "(", strand(query_blast_subset_ranges), ")")
  
  # combine sequences with original query sequence
  query_blast_subset_seq <- c(query_seq, query_blast_subset_seq)
  
  # write compiled sequences to temporary fasta file
  Biostrings::writeXStringSet(x = query_blast_subset_seq, filepath = paste0(out_dir, "temp.fa"), append = F)
  query_repeat_subbed <- gsub(".fasta", "", repeat_name) 
  # create alignment of sequences, file named after original query with "/" character replaced with "_"
  system(paste0("mafft --localpair --thread 12 ", out_dir, "temp.fa > ", out_dir, species_name, "_", query_repeat_subbed, "_aligned.fa"))

