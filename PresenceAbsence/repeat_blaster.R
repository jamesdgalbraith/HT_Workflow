# Script to search for repeats in other metazoan genomes (relaxed blast parameters)

suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

repeats <- c("~/HT_Workflow/Proto2-Snek/Proto2-Snek.fasta", "~/HT_Workflow/RTE-Snek/RTE-Snek.fasta", "~/HT_Workflow/Rex1-Snek_1/Rex1-Snek_1.fasta", "~/HT_Workflow/Rex1-Snek_2/Rex1-Snek_2.fasta")


species_list <- read_tsv("~/all_genomes.txt", col_names = c("clade_name", "species_name", "genome_name")) %>%
  mutate(genome_path = paste0("~/Genomes/", clade_name, "/", species_name, "/", genome_name))

for(j in 1:length(repeats)){
  
  query_repeat <- repeats[j]
  
  contains_TE <- tibble(genome_name = character())
  
  all_blast_out <- tibble(seqnames = character(), sstart = integer(), send = integer(), pident = double(), qcovs = integer(), bitscore = double(), length = integer(), species_name = character())
  
  for(i in 1:nrow(species_list)){
    
    genome_path <- species_list$genome_path[i]
    
    print(paste0("Searching for ", query_repeat, " in ", species_list$genome_path[i]))
    
    blast_out <- read.table(text=system(paste0("blastn -evalue 0.00002 -num_threads 12 -reward 3 -penalty -4 -xdrop_ungap 80 -xdrop_gap 130 -xdrop_gap_final 150 -word_size 7 -dust yes -gapopen 30 -gapextend 6 -query ", query_repeat, " -db ", genome_path, " -outfmt \"6 sseqid sstart send pident qcovs bitscore length\""), intern = TRUE), col.names = c("seqnames", "sstart", "send", "pident", "qcovs", "bitscore", "length")) %>%
      as_tibble() %>%
      filter(length > 500, pident > 70) %>%
      mutate(seqnames = as.character(seqnames))
    
    if(nrow(blast_out) > 0){
      contains_TE <- rbind(contains_TE, tibble(genome_name = genome_path))
      
      all_blast_out <- blast_out %>%
        mutate(species_name = genome_path) %>%
        rbind(all_blast_out)
    }
    
  }
  
  write_tsv(all_blast_out, paste0("blaster/", sub("~.*/", "", sub(".fasta", "", query_repeat)), "_in_other_species.tsv"), col_names = T)
  
}

# Count number in each species
Proto2 <- read_tsv("blaster/Proto2-Snek_in_other_species.tsv", col_names = c("sseqid", "sstart", "send", "pident", "qcovs", "bitscore", "length", "genome_path"), skip = 1) %>%
  select(genome_path) %>%
  table()
Proto2 <- tibble(genome = names(Proto2), Proto2 = as.integer(Proto2))

Rex_1 <- read_tsv("blaster/Rex1-Snek_1_in_other_species.tsv", col_names = c("sseqid", "sstart", "send", "pident", "qcovs", "bitscore", "length", "genome_path"), skip = 1) %>%
  select(genome_path) %>%
  table()
Rex_1 <- tibble(genome = names(Rex_1), Rex_1 = as.integer(Rex_1))

Rex_2 <- read_tsv("blaster/Rex1-Snek_2_in_other_species.tsv", col_names = c("sseqid", "sstart", "send", "pident", "qcovs", "bitscore", "length", "genome_path"), skip = 1) %>%
  select(genome_path) %>%
  table()
Rex_2 <- tibble(genome = names(Rex_2), Rex_2 = as.integer(Rex_2))

RTE <- read_tsv("blaster/RTE-Snek_in_other_species.tsv", col_names = c("sseqid", "sstart", "send", "pident", "qcovs", "bitscore", "length", "genome_path"), skip = 1) %>%
  select(genome_path) %>%
  table()
RTE <- tibble(genome = names(RTE), RTE = as.integer(RTE))


all_joined <- full_join(Proto2, Rex_1) %>%
  full_join(Rex_2) %>%
  full_join(RTE)

all_joined <- all_joined %>%
  replace_na(list(Proto2 = 0, Rex_1 = 0, Rex_2 = 0, RTE = 0))

all_joined <- all_joined %>%
  separate(genome, into = c("tilde", "genomes", "clade", "Species", "fasta"), sep = "/") %>%
  select(Species, Proto2, Rex_1, Rex_2, RTE) %>%
  dplyr::rename(`Proto2-Snek` = Proto2, `Rex1-Snek_1` = Rex_1, `Rex1-Snek_2` = Rex_2, `RTE-Snek` = RTE)

write_tsv(x = all_joined, path = "species_repeat_count.tsv")
