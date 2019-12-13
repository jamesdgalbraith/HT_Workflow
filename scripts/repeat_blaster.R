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
j=1
query_repeat <- repeats[j]

blast_out <- read_tsv(paste0("blaster/word_size_10/", sub("~.*/", "", sub(".fasta", "", query_repeat)), "_in_other_species.tsv"), col_names = T) %>%
  mutate(species_name = gsub("/", ",", species_name, )) %>%
  separate(species_name, into = c("blank", "home", "james", "Genomes", "clade_name", "species_name", "genome_name"), sep = ",") %>%
  dplyr::rename(start = sstart, end = send) %>%
  select(-blank, -home, -james, -Genomes)

species_hits <- blast_out %>%
  select(clade_name, species_name, genome_name) %>%
  base::unique() %>%
  arrange(species_name)

blast_out <- blast_out %>%
  arrange(-pident)

