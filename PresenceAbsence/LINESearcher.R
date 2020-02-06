#!/usr/bin/env Rscript

## Script  for creating initial alignements for manual annotation of repeats

suppressMessages(library(tidyverse))

species_hits <- read_tsv("~/all_genomes.txt", col_names = c("clade_name", "species_name", "genome_name"))
repeat_names <- c("Rex1-Snek_1", "Rex1-Snek_2", "Rex1-Snek_3", "Rex1-Snek_4", "RTE-Snek", "Proto2-Snek")


# set output folder

for(j in 1:base::length(repeat_names)){
  
  query_repeat <- paste0("~/HT_Workflow/", repeat_names[j], "/", repeat_names[j], ".fasta")
  
  for(i in 1:nrow(species_hits)){
    clade <- species_hits$clade_name[i]
    species_name <- species_hits$species_name[i]
    genome_name <- species_hits$genome_name[i]
    print(species_name)
    
    # set and read in genome and index
    genome_path <- paste0("~/Genomes/", clade, "/", species_name, "/", genome_name)
    
    blast_out <- read.table(text=system(paste0("blastn -evalue 0.00002 -num_threads 12 -word_size 10 -dust yes -query ", query_repeat, " -db ", genome_path, " -outfmt 6 -reward 3 -penalty -4 -xdrop_ungap 80 -xdrop_gap 130 -xdrop_gap_final 150 -gapopen 30 -gapextend 6"), intern = TRUE), col.names = c("sseqid", "qseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend","sstart", "send", "evalue", "bitscore")) %>% as_tibble()
    
    blast_out <- blast_out %>%
      filter(length > 1000)
    
    # if no reuslts skip to next species
    if(nrow(blast_out) < 1){
      next
    } else {
      write_tsv(blast_out, paste0("~/HT_Workflow/new_method/", repeat_names[j], "/", species_name, ".tsv"), col_names = F)
    }
    
  }
}
