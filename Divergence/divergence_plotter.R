# Script for plotting Jukes Cantor divergence of repeats found in sea snakes

library(tidyverse)

repeat_list <- tibble(repeat_name = c("Proto2-Snek", "Rex1-Snek_1", "Rex1-Snek_2", "Rex1-Snek_3", "Rex1-Snek_4", "RTE-Snek", "RTE-Snek_2"), max = c(260, 50, 125, 200, 80, 110, 150))
species_list <- tibble(species = c("Aipysurus_laevis", "Emydocephalus_ijimae", "Hydrophis_melanocephalus", "Notechis_scutatus", "Pseudonaja_textilis", "Laticauda_colubrina", "Ophiophagus_hannah"),
                       genome = c("assembly_20171114.fasta", "emyIji_1.0.fasta", "hydMel_1.0.fasta", "TS10Xv2-PRI.fasta", "EBS10Xv2-PRI.fasta", "latCor_1.0.fasta", "OphHan1.0.fasta"))

# i=5
# j=2
for(i in 1:nrow(species_list)){
  for(j in 1:nrow(repeat_list)){
    # set query repeat and genome path
    query_repeat <- paste0("~/HT_Workflow/", "Divergence", "/", repeat_list$repeat_name[j], ".fasta")
    genome_path <- paste0("~/Genomes/Reptiles/", species_list$species[i], "/", species_list$genome[i])
    
    # perform blastn search
    blastn <- read.table(text=system(paste0("blastn -dust yes -query ", query_repeat, " -db ", genome_path, " -outfmt \"6 sseqid sstart send pident qcovs bitscore length mismatch evalue\""), intern = TRUE), col.names = c("seqnames", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue")) %>%
      as_tibble() %>%
      dplyr::filter(length >= 50, pident > 94) %>%
      mutate(div = 100 - pident, d = mismatch/length, jc_dist = (-3 / 4) * log(1 - (4 * d / 3))) # calculate Jukes-Cantor distance
    
    # plot resulting graphs (max y axes determined manually)
    ggplot(blastn, aes(d)) + geom_histogram(binwidth = 0.01) +
      scale_x_continuous(limits = c(-0.01,0.26), expand = c(0,0)) + labs(x = NULL, y = NULL) + scale_y_continuous(limits = c(0, repeat_list$max[j]), expand = c(0,0))
    
    # save plots
    ggsave(filename = paste0("~/HT_Workflow/Divergence/JCdist/", repeat_list$repeat_name[j], "_in_", species_list$species[i], "_JCdist.pdf"), width = 10, height = 10, units = "cm")
  }
}
warnings()

for(i in 1:nrow(species_list)){
  for(j in 1:nrow(repeat_list)){
    # set query repeat and genome path
    query_repeat <- paste0("~/HT_Workflow/", "Divergence", "/", repeat_list$repeat_name[j], ".fasta")
    genome_path <- paste0("~/Genomes/Reptiles/", species_list$species[i], "/", species_list$genome[i])
    
    # perform blastn search
    blastn <- read.table(text=system(paste0("blastn -dust yes -query ", query_repeat, " -db ", genome_path, " -outfmt \"6 sseqid sstart send pident qcovs bitscore length mismatch evalue\""), intern = TRUE), col.names = c("seqnames", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue")) %>%
      as_tibble() %>%
      dplyr::filter(length >= 50) %>%
      mutate(div = 100 - pident, d = mismatch/length, jc_dist = (-3 / 4) * log(1 - (4 * d / 3))) # calculate Jukes-Cantor distance
    
    # plot resulting graphs (max y axes determined manually)
    ggplot(blastn, aes(div)) + geom_histogram(binwidth = 1) +
      scale_x_continuous(limits = c(-1, 30), expand = c(0,0)) + labs(x = NULL, y = NULL) + scale_y_continuous(limits = c(0, repeat_list$max[j]), expand = c(0,0))
    
    # save plots
    ggsave(filename = paste0("~/HT_Workflow/Divergence/JCdist/", repeat_list$repeat_name[j], "_in_", species_list$species[i], "_pident.pdf"), width = 10, height = 10, units = "cm")
  }
}
warnings()
