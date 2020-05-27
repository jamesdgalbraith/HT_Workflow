# Script for plotting Jukes Cantor divergence of repeats found in sea snakes

library(BSgenome)
library(GenomicRanges)
library(plyranges)
library(tidyverse)

if(!dir.exists(paste0("~/HT_Workflow/HT_Workflow/Divergence/scratch/"))){
  dir.create(paste0("~/HT_Workflow/HT_Workflow/Divergence/scratch/"), recursive = T)
  }

repeat_list <- tibble(repeat_name = c("Proto2-Snek", "Rex1-Snek_1", "Rex1-Snek_2", "Rex1-Snek_3", "Rex1-Snek_4", "RTE-Snek_1", "RTE-Snek_2", "RTE-Kret"),
                      count_max = c(300, 30, 120, 80, 100, 125, 400, 120),
                      bp_max = c(120000, 10000, 20000, 17000, 30000, 41000, 150000, 200000))

species_list <- tibble(species = c("Aipysurus_laevis", "Emydocephalus_ijimae", "Hydrophis_melanocephalus", "Notechis_scutatus", "Pseudonaja_textilis", "Laticauda_colubrina", "Ophiophagus_hannah"),
                       genome = c("assembly_20171114.fasta", "emyIji_1.0.fasta", "hydMel_1.0.fasta", "TS10Xv2-PRI.fasta", "EBS10Xv2-PRI.fasta", "latCor_1.0.fasta", "OphHan1.0.fasta"))

for(i in 1:nrow(species_list)){
  genome_path <- paste0("~/Genomes/Reptiles/", species_list$species[i], "/", species_list$genome[i])
  
  blastn <- read.table(text=system(paste0("blastn -dust yes -query ~/HT_Workflow/HT_Workflow/compiled_LINEs.fasta -db ", genome_path, " -outfmt \"6 sseqid sstart send pident qcovs bitscore length mismatch evalue qseqid\""), intern = TRUE), col.names = c("seqnames", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "qseqid")) %>%
    as_tibble() %>%
    dplyr::filter(length >= 50) %>%
    mutate(pident_rnd = ceiling(pident), div = 100 - pident, d = mismatch/length, seqnames = as.character(seqnames), qseqid = as.character(qseqid),
           jc_dist = (-3 / 4) * log(1 - (4 * d / 3))) # calculate Jukes-Cantor distance
  
  # read in genome and garbage collect
  genome_seq <- Biostrings::readDNAStringSet(filepath = genome_path)
  gc()
  names(genome_seq) <- sub(" .*", "", names(genome_seq))
  
  
  for(j in 1:nrow(repeat_list)){
    
    query_repeat <- paste0("~/HT_Workflow/", repeat_list$repeat_name[j], "/", repeat_list$repeat_name[j], ".fasta")
    
    bed <- blastn %>%
      filter(qseqid == repeat_list$repeat_name[j]) %>%
      mutate(strand = case_when(sstart > send ~ "-", send > sstart ~ "+"),
             start = case_when(sstart < send ~ sstart, send < sstart ~ send),
             end = case_when(sstart > send ~ sstart, send > sstart ~ send)) %>%
      select(seqnames, start, end, strand) %>%
      as_granges() %>%
      GenomicRanges::reduce(min.gapwidth = 200)
    
    if(nrow(as_tibble(bed))>1){
      
      # get seqs
      seqs <- Biostrings::getSeq(genome_seq, bed)
      
      names(seqs) <- paste0(seqnames(bed), ":", ranges(bed), "(", strand(bed), ")")

      Biostrings::writeXStringSet(x = seqs, filepath = paste0("~/HT_Workflow/HT_Workflow/Divergence/scratch/temp.fa"))
      
      blastn_recip <- read_tsv(system(paste0("blastn -dust yes -query ~/HT_Workflow/HT_Workflow/Divergence/scratch/temp.fa -subject ~/HT_Workflow/HT_Workflow/compiled_LINEs.fasta -outfmt \"6 qseqid sseqid qstart qend sstart send pident qcovs bitscore length evalue qlen slen\""), intern = TRUE), col_names = c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "pident", "qcovs", "bitscore", "length", "evalue", "qlen", "slen")) %>%
        dplyr::filter(length >= 50) %>%
        mutate(pident_rnd = ceiling(pident), div = 100 - pident_rnd)
      
      blastn_recip <- blastn_recip %>%
        group_by(qseqid) %>%
        arrange(-pident_rnd) %>%
        dplyr::slice(1) %>%
        ungroup() %>%
        filter(sseqid == as.character(repeat_list$repeat_name[j]))
      
      blastn_recip %>%
        filter(pident >= 97, length >= 50)
      
      
      blastn_recip_bp <- blastn_recip %>%
        group_by(div) %>%
        select(div, length) %>%
        summarise(bp = sum(length))
    } else {
      blastn_recip <- tibble(div = double(), bp = double())
      blastn_recip_bp <- tibble(div = 0:30, bp = 0)
    }
    
    
    # create plot of insertions vs divergence
    ggplot(blastn_recip, aes(div)) + geom_bar(binwidth = 1) + scale_x_continuous(expand = c(0,0), limits = c(-1, 31))  + labs(x = NULL, y = NULL) +
    scale_y_continuous(limits = c(0, repeat_list$count_max[j]), expand = c(0,0)) + theme_bw()
    ggsave(filename = paste0("~/HT_Workflow/Divergence/plots/", repeat_list$repeat_name[j], "_in_", species_list$species[i], "_counts.pdf"), width = 10, height = 10, units = "cm")

    # create plot of base pairs vs divergence
    ggplot(blastn_recip_bp, aes(div, bp)) + geom_bar(stat = "identity", width=1, position=position_dodge(width=0)) + scale_x_continuous(expand = c(0,0), limits = c(-1, 31)) +
      labs(x = NULL, y = NULL) + scale_y_continuous(limits = c(0, repeat_list$bp_max[j]), expand = c(0,0)) + theme_bw()
    ggsave(filename = paste0("~/HT_Workflow/Divergence/plots/", repeat_list$repeat_name[j], "_in_", species_list$species[i], "_bp.pdf"), width = 10, height = 10, units = "cm")
    
    # filter out hits below 97% id
    blastn_recip <- blastn_recip %>%
    filter(pident_rnd >= 97)


    if(nrow(blastn_recip)>0){
    # select full length sequences
    blastn_recip_fl <- blastn_recip %>%
      filter(length >= slen * 0.9)

    # set length of consensus repeat
    cons_len <- blastn_recip$slen[1]

    # create matrix for coverage
    bp_coverage <- matrix(rep(0, length(blastn_recip$sseqid)*as.numeric(cons_len)), byrow = T, ncol = as.numeric(cons_len))

    # calculate coverage of each base pair
    for(k in 1:length(blastn_recip$sseqid)){
      bp_coverage[k,]<-c(rep(0,blastn_recip$sstart[k]-1),rep(1,abs(blastn_recip$send[k]-blastn_recip$sstart[k])+1), rep(0,as.numeric(cons_len)-blastn_recip$send[k]))
    }

    # convert matrix to tibble
    bp_coverage_tbl <- tibble(bp_x = 1:cons_len, bp_y = colSums(bp_coverage)) %>%
      mutate(bp_y_scaled = 3 * bp_y / max(bp_y))

    # plot coverage
    ggplot() +
      geom_segment(data = blastn_recip, aes(x = sstart, xend = send, y = 100 - pident, yend = 100 - pident), colour = "black") +
      geom_segment(data = blastn_recip_fl, aes(x = sstart, xend = send, y = 100 - pident, yend = 100 - pident), colour = "#db7800", size = 1.5) +
      geom_line(data = bp_coverage_tbl, aes(bp_x, bp_y_scaled), colour = "#006edb") +
      scale_y_continuous(limits = c(-0.05,3.05), sec.axis = sec_axis(trans = ~ . * max(bp_coverage_tbl$bp_y)/3), expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      theme_bw() + labs(x = NULL, y = NULL)

    ggsave(filename = paste0("~/HT_Workflow/Divergence/plots/", repeat_list$repeat_name[j], "_in_", species_list$species[i], "_coverage.pdf"), width = 10, height = 10, units = "cm")}

  }
}

    
# genome_seq <- NULL
gc()

