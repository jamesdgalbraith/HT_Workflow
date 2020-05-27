#!/usr/bin/env Rscript

## Script  for creating initial alignements for manual annotation of repeats

suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(BSgenome))

flank5 <- 0
flank3 <- 0

species_hits <- read_tsv("~/all_genomes.txt", col_names = c("clade_name", "species_name", "genome_name"))

repeat_names <- c("Proto2-Snek", "Rex1-Snek_1", "Rex1-Snek_2", "Rex1-Snek_3", "Rex1-Snek_4", "RTE-Snek_1", "RTE-Snek_2")

# set output folder
for(j in 1:base::length(repeat_names)){
  
  query_repeat <- paste0("~/HT_Workflow/", repeat_names[j], "/", repeat_names[j], ".fasta")
  repeat_seq <- readDNAStringSet(query_repeat)
  
  if(file.exists(paste0("PresenceAbsence/curation/", repeat_names[j], "/species_info.tsv"))){
    file.remove(paste0("PresenceAbsence/curation/", repeat_names[j], "/species_info.tsv"))
    }
  
  out_folder <- paste0("PresenceAbsence/curation/", repeat_names[j], "/")
  
  if(!dir.exists(out_folder)){dir.create(out_folder, recursive = T)}
  
  for(i in 1:nrow(species_hits)){
    clade <- species_hits$clade_name[i]
    species_name <- species_hits$species_name[i]
    genome_name <- species_hits$genome_name[i]
    print(species_name)
    
    # set and read in genome and index
    genome_path <- paste0("~/Genomes/", clade, "/", species_name, "/", genome_name)
    
    blast_out <- read.table(text=system(paste0("blastn -evalue 0.00002 -num_threads 12 -word_size 10 -dust yes -query ", query_repeat, " -db ", genome_path, " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\" -reward 3 -penalty -4 -xdrop_ungap 80 -xdrop_gap 130 -xdrop_gap_final 150 -gapopen 30 -gapextend 6"), intern = TRUE), col.names = c("qseqid", "seqnames", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen")) %>% as_tibble()
    
    blast_out <- blast_out %>%
      filter(length > 400)
    
    # if no reuslts skip to next species
    if(nrow(blast_out) < 1){
      next
    } else {
      # write output to file
      
      write_tsv(x = species_hits[i,], path = paste0("PresenceAbsence/curation/", repeat_names[j], "/species_info.tsv"), append = T, col_names = F)
      
      # import genome
      genome_path <- paste0("~/Genomes/", clade, "/", species_name, "/", genome_name, "")
      genome_seq <- Biostrings::readDNAStringSet(filepath = genome_path)
      gc()
      names(genome_seq) <- sub(" .*", "", names(genome_seq))
      genome_fai <- tibble(seqnames = names(genome_seq), scaffold_length = as.double(width(genome_seq)))
      
      bed <- blast_out %>%
        as_tibble() %>%
        dplyr::mutate(seqnames = as.character(seqnames), sstart = as.double(sstart), send = as.double(send), length = as.double(length)) %>%
        dplyr::arrange(-bitscore) %>%
        mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+"),
               start = case_when(sstart < send ~ sstart, send < sstart ~ send),
               end = case_when(sstart > send ~ sstart, send > sstart ~ send)) %>%
        inner_join(genome_fai) %>%
        arrange(-bitscore) %>%
        mutate(start = case_when(strand == "+" ~ start - flank5, strand == "-" ~ start - flank3),
               end = case_when(strand == "-" ~ end + flank5, strand == "+" ~ end + flank3)) %>%
        mutate(start = case_when(start <= 1 ~ 1, start > 1 ~ start),
               end = case_when(end > scaffold_length ~ scaffold_length, end <= scaffold_length ~ end)) %>%
        mutate(name = paste0(seqnames, ":", start, "-", end, "(", strand, ")")) %>%
        dplyr::select(seqnames, start, end, strand, name)
      
      if(nrow(bed)>25){
        bed <- bed %>% 
          dplyr::slice(1:25)
      }
      
      bed <- plyranges::as_granges(bed)  
    
      seqs <- Biostrings::getSeq(genome_seq, bed)
      
      # name seqs
      names(seqs) <- bed$name
      
      # merge with original query
      seqs <- c(repeat_seq, seqs)
      
      # write seqs to files
      Biostrings::writeXStringSet(x = seqs, filepath = paste0(out_folder, species_name, "_", repeat_names[j], "_hits.fa"))
      
      # perform multiple alignment
      system(paste0("mafft --localpair --maxiterate 10 --thread 12 ", out_folder, species_name, "_", repeat_names[j], "_hits.fa > ", out_folder, species_name, "_", repeat_names[j], "_aligned.fa"))
        
    }
    
  }
}
