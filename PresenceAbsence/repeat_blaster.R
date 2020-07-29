# Script to search for repeats in other metazoan genomes (relaxed blast parameters)

suppressMessages(library(tidyverse))

# set repeat names
repeats <- c("Proto2-Snek", "RTE-Snek_1", "RTE-Snek_2", "Rex1-Snek_1", "Rex1-Snek_2", "Rex1-Snek_3", "Rex1-Snek_4", "RTE-Kret")

species_list <- read_tsv("~/all_genomes.txt", col_names = c("clade_name", "species_name", "genome_name")) %>%
  mutate(genome_path = paste0("~/Genomes/", clade_name, "/", species_name, "/", genome_name))

for(j in 1:length(repeats)){
  
  query_repeat <- repeats[j]
  
  contains_TE <- tibble(genome_name = character())
  
  all_blast_out <- tibble(seqnames = character(), sstart = integer(), send = integer(), pident = double(), qcovs = integer(), bitscore = double(), length = integer(), species_name = character())
  
  for(i in 1:nrow(species_list)){
    
    genome_path <- species_list$genome_path[i]
    
    print(paste0("Searching for ", query_repeat, " in ", species_list$genome_path[i]))
    
    blast_out <- read.table(text=system(paste0("blastn -evalue 0.00002 -num_threads 12 -reward 3 -penalty -4 -xdrop_ungap 80 -xdrop_gap 130 -xdrop_gap_final 150 -word_size 7 -dust yes -gapopen 30 -gapextend 6 -query ~/HT_Workflow/", query_repeat, "/", query_repeat, ".fasta -db ", genome_path, " -outfmt \"6 sseqid sstart send pident qcovs bitscore length\""), intern = TRUE), col.names = c("seqnames", "sstart", "send", "pident", "qcovs", "bitscore", "length")) %>%
      as_tibble() %>%
      filter(length > 500) %>%
      mutate(seqnames = as.character(seqnames))
    
    if(nrow(blast_out) > 0){
      contains_TE <- rbind(contains_TE, tibble(genome_name = genome_path))
      
      all_blast_out <- blast_out %>%
        mutate(species_name = genome_path) %>%
        rbind(all_blast_out)
    }
    
  }
  
  write_tsv(all_blast_out, paste0("PresenceAbsence/blaster/",  query_repeat, "_in_other_species.tsv"), col_names = T)
  
}

# Count number in each species
Proto2_out <- read_tsv("PresenceAbsence/blaster/Proto2-Snek_in_other_species.tsv", col_names = c("sseqid", "sstart", "send", "pident", "qcovs", "bitscore", "length", "genome_path"), skip = 1) %>%
  filter(length >= 1000)

Proto2_above_75 <- Proto2_out %>%
  filter(pident >= 75) %>%
  select(genome_path)

Proto2_above_75 <- tibble(species_name = names(table(Proto2_above_75)), Proto2 = as.integer(table(Proto2_above_75)))

Proto2_below_75 <- Proto2_out %>%
  filter(!genome_path %in% Proto2_above_75$species_name) %>%
  select(genome_path)

Proto2_below_75 <- tibble(species_name = names(table(Proto2_below_75)), Proto2 = as.integer(table(Proto2_below_75)))

RTE_1_out <- read_tsv("PresenceAbsence/blaster/RTE-Snek_1_in_other_species.tsv", col_names = c("sseqid", "sstart", "send", "pident", "qcovs", "bitscore", "length", "genome_path"), skip = 1) %>%
  filter(length >= 1000)

RTE_1_above_75 <- RTE_1_out %>%
  filter(pident >= 75) %>%
  select(genome_path)

RTE_1_above_75 <- tibble(species_name = names(table(RTE_1_above_75)), RTE_1 = as.integer(table(RTE_1_above_75)))

RTE_1_below_75 <- RTE_1_out %>%
  filter(!genome_path %in% RTE_1_above_75$species_name) %>%
  select(genome_path)

RTE_1_below_75 <- tibble(species_name = names(table(RTE_1_below_75)), RTE_1 = as.integer(table(RTE_1_below_75)))

RTE_2_out <- read_tsv("PresenceAbsence/blaster/RTE-Snek_2_in_other_species.tsv", col_names = c("sseqid", "sstart", "send", "pident", "qcovs", "bitscore", "length", "genome_path"), skip = 1) %>%
  filter(length >= 1000)

RTE_2_above_75 <- RTE_2_out %>%
  filter(pident >= 75) %>%
  select(genome_path)

RTE_2_above_75 <- tibble(species_name = names(table(RTE_2_above_75)), RTE_2 = as.integer(table(RTE_2_above_75)))

RTE_2_below_75 <- RTE_2_out %>%
  filter(!genome_path %in% RTE_2_above_75$species_name) %>%
  select(genome_path)

RTE_2_below_75 <- tibble(species_name = names(table(RTE_2_below_75)), RTE_2 = as.integer(table(RTE_2_below_75)))

RTE_Kret_out <- read_tsv("PresenceAbsence/blaster/RTE-Kret_in_other_species.tsv", col_names = c("sseqid", "sstart", "send", "pident", "qcovs", "bitscore", "length", "genome_path"), skip = 1) %>%
  filter(length >= 1000)

RTE_Kret_above_75 <- RTE_Kret_out %>%
  filter(pident >= 75) %>%
  select(genome_path)

RTE_Kret_above_75 <- tibble(species_name = names(table(RTE_Kret_above_75)), RTE_Kret = as.integer(table(RTE_Kret_above_75)))

RTE_Kret_below_75 <- RTE_Kret_out %>%
  filter(!genome_path %in% RTE_Kret_above_75$species_name) %>%
  select(genome_path)

RTE_Kret_below_75 <- tibble(species_name = names(table(RTE_Kret_below_75)), RTE_Kret = as.integer(table(RTE_Kret_below_75)))

Rex1_1_out <- read_tsv("PresenceAbsence/blaster/Rex1-Snek_1_in_other_species.tsv", col_names = c("sseqid", "sstart", "send", "pident", "qcovs", "bitscore", "length", "genome_path"), skip = 1) %>%
  filter(length >= 1000)

Rex1_1_above_75 <- Rex1_1_out %>%
  filter(pident >= 75) %>%
  select(genome_path)

Rex1_1_above_75 <- tibble(species_name = names(table(Rex1_1_above_75)), Rex1_1 = as.integer(table(Rex1_1_above_75)))

Rex1_1_below_75 <- Rex1_1_out %>%
  filter(!genome_path %in% Rex1_1_above_75$species_name) %>%
  select(genome_path) %>%
  base::unique()

Rex1_1_below_75 <- tibble(species_name = names(table(Rex1_1_below_75)), Rex1_1 = as.integer(table(Rex1_1_below_75)))

Rex1_2_out <- read_tsv("PresenceAbsence/blaster/Rex1-Snek_2_in_other_species.tsv", col_names = c("sseqid", "sstart", "send", "pident", "qcovs", "bitscore", "length", "genome_path"), skip = 1) %>%
  filter(length >= 1000)

Rex1_2_above_75 <- Rex1_2_out %>%
  filter(pident >= 75) %>%
  select(genome_path)

Rex1_2_above_75 <- tibble(species_name = names(table(Rex1_2_above_75)), Rex1_2 = as.integer(table(Rex1_2_above_75)))

Rex1_2_below_75 <- Rex1_2_out %>%
  filter(!genome_path %in% Rex1_2_above_75$species_name) %>%
  select(genome_path)

Rex1_2_below_75 <- tibble(species_name = names(table(Rex1_2_below_75)), Rex1_2 = as.integer(table(Rex1_2_below_75)))

Rex1_3_out <- read_tsv("PresenceAbsence/blaster/Rex1-Snek_3_in_other_species.tsv", col_names = c("sseqid", "sstart", "send", "pident", "qcovs", "bitscore", "length", "genome_path"), skip = 1) %>%
  filter(length >= 1000)

Rex1_3_above_75 <- Rex1_3_out %>%
  filter(pident >= 75) %>%
  select(genome_path)

Rex1_3_above_75 <- tibble(species_name = names(table(Rex1_3_above_75)), Rex1_3 = as.integer(table(Rex1_3_above_75)))

Rex1_3_below_75 <- Rex1_3_out %>%
  filter(!genome_path %in% Rex1_3_above_75$species_name) %>%
  select(genome_path)

Rex1_3_below_75 <- tibble(species_name = names(table(Rex1_3_below_75)), Rex1_3 = as.integer(table(Rex1_3_below_75)))

Rex1_4_out <- read_tsv("PresenceAbsence/blaster/Rex1-Snek_4_in_other_species.tsv", col_names = c("sseqid", "sstart", "send", "pident", "qcovs", "bitscore", "length", "genome_path"), skip = 1) %>%
  filter(length >= 1000)

Rex1_4_above_75 <- Rex1_4_out %>%
  filter(pident >= 75) %>%
  select(genome_path)

Rex1_4_above_75 <- tibble(species_name = names(table(Rex1_4_above_75)), Rex1_4 = as.integer(table(Rex1_4_above_75)))

Rex1_4_below_75 <- Rex1_4_out %>%
  filter(!genome_path %in% Rex1_4_above_75$species_name) %>%
  select(genome_path)

Rex1_4_below_75 <- tibble(species_name = names(table(Rex1_4_below_75)), Rex1_4 = as.integer(table(Rex1_4_below_75)))

table_1 <- Rex1_1_above_75 %>%
  full_join(Rex1_2_above_75) %>%
  full_join(Rex1_3_above_75) %>%
  full_join(Rex1_4_above_75) %>%
  full_join(RTE_1_above_75) %>%
  full_join(RTE_2_above_75) %>%
  full_join(RTE_Kret_above_75) %>%
  full_join(Proto2_above_75)
  
table_2 <- table_1 %>%
  arrange(species_name) %>%
  mutate(RTE_1 = ifelse(is.na(RTE_1), 0, RTE_1),
         RTE_2 = ifelse(is.na(RTE_2), 0, RTE_2),
         RTE_Kret = ifelse(is.na(RTE_Kret), 0, RTE_Kret),
         Proto2 = ifelse(is.na(Proto2), 0, Proto2),
         Rex1_1 = ifelse(is.na(Rex1_1), 0, Rex1_1),
         Rex1_2 = ifelse(is.na(Rex1_2), 0, Rex1_2),
         Rex1_3 = ifelse(is.na(Rex1_3), 0, Rex1_3),
         Rex1_4 = ifelse(is.na(Rex1_4), 0, Rex1_4),
         species_name = sub("~/Genomes/", "", species_name)) %>%
  separate(species_name, into = c("clade", "species_name", "genome_name"), sep = "/") %>%
  select(-genome_name)

write_tsv(x = table_2, path = "PresenceAbsence/blaster/above_75.tsv")

table_3 <- Rex1_1_below_75 %>%
  full_join(Rex1_2_below_75) %>%
  full_join(Rex1_3_below_75) %>%
  full_join(Rex1_4_below_75) %>%
  full_join(RTE_1_below_75) %>%
  full_join(RTE_2_below_75) %>%
  full_join(RTE_Kret_below_75) %>%
  full_join(Proto2_below_75)

table_4 <- table_3 %>%
  arrange(species_name) %>%
  mutate(RTE_1 = ifelse(is.na(RTE_1), 0, RTE_1),
         RTE_2 = ifelse(is.na(RTE_2), 0, RTE_2),
         RTE_Kret = ifelse(is.na(RTE_Kret), 0, RTE_Kret),
         Proto2 = ifelse(is.na(Proto2), 0, Proto2),
         Rex1_1 = ifelse(is.na(Rex1_1), 0, Rex1_1),
         Rex1_2 = ifelse(is.na(Rex1_2), 0, Rex1_2),
         Rex1_3 = ifelse(is.na(Rex1_3), 0, Rex1_3),
         Rex1_4 = ifelse(is.na(Rex1_4), 0, Rex1_4),
         species_name = sub("~/Genomes/", "", species_name)) %>%
  separate(species_name, into = c("clade", "species_name", "genome_name"), sep = "/") %>%
  select(-genome_name)

write_tsv(x = table_4, path = "PresenceAbsence/blaster/below_75.tsv")