library(tidyverse)
install.packages("svglite")

# set query and genome
genome_path <- "~/Genomes/Reptiles/Laticauda_colubrina/latCor_1.0.fasta"
query_repeat <- "~/HT_Workflow/RTE-Kret/RTE-Kret.fasta"

# perform blastn search
blastn <- read_tsv(system(paste0("blastn -num_threads 4 -dust yes -query ", query_repeat, " -db ", genome_path, " -outfmt \"6 sseqid qseqid qstart qend sstart send pident qcovs bitscore length mismatch evalue qlen slen\" -task megablast"), intern = TRUE), col_names = c("seqnames", "qseqid",  "qstart", "qend", "sstart", "send", "pident", "qcovs", "bitscore", "length", "mismatch", "evalue", "qlen", "slen")) %>%
  as_tibble() %>%
  dplyr::filter(length >= 50) %>%
  mutate(div = 100 - pident, d = mismatch/length, jc_dist = (-3 / 4) * log(1 - (4 * d / 3)), jc_dist_round = as.integer(100 * round(jc_dist, digits = 2))) %>%
  arrange(-bitscore)

query_width <- blastn$qlen[1]

# select near full length repeats
blastn_fl <- blastn %>%
  filter(length >= qlen * 0.9)

# calculate coverage of each bp
bp_coverage <- matrix(rep(0, length(blastn$seqnames)*as.numeric(query_width)), byrow = T, ncol = as.numeric(query_width))

for(k in 1:length(blastn$seqnames)){
  bp_coverage[k,]<-c(rep(0,blastn$qstart[k]-1),rep(1,abs(blastn$qend[k]-blastn$qstart[k])+1), rep(0,as.numeric(query_width)-blastn$qend[k]))
}

# convert matrix to tibble
bp_coverage_tbl <- tibble(bp_x = 1:query_width, bp_y = colSums(bp_coverage)) %>%
  mutate(bp_y_scaled = max(blastn$div) * bp_y / max(bp_y))

blastn %>%
  filter(qend < 2080)

# plot coverage
ggplot() + geom_segment(data = blastn, aes(x = qstart, xend = qend, y = 100 - pident, yend = 100 - pident), colour = "black") + # plot all repeats
  
  geom_segment(data = blastn_fl, aes(x = qstart, xend = qend, y = 100 - pident, yend = 100 - pident), colour = "#db7800", size = 1.5) + # plot full length repeats
  geom_line(data = bp_coverage_tbl, aes(bp_x, bp_y_scaled), colour = "#006edb", size = 1.5) + # plot coverage data
  scale_y_continuous(sec.axis = sec_axis(trans = ~ . * max(bp_coverage_tbl$bp_y)/3, name = "coverage"), # add second axis
                     limits = c(-(max(blastn$div) *0.01),(max(blastn$div) *1.01)),  # scale for second axis
                     expand = c(0,0)) + # don't expand axis
  scale_x_continuous(expand = c(0,0)) + # don't expand axis
  theme_bw() +
  labs(x = "bp", y = "pairwise divergence from consensus", title = blastn$qseqid[1])

ggsave("TreeBuilding/RTE-Kret_coverage.svg", device = "svg", height = 10, width = 10, units = "cm")

# plot coverage
jc_dist_mat <- matrix(nrow = length(blastn$seqnames), ncol = max(blastn$jc_dist_round) + 1)

for(k in 1:length(blastn$seqnames)){
  jc_dist_mat[k,] <-c(rep(0,blastn$jc_dist_round[k]),blastn$length[k], rep(0,max(blastn$jc_dist_round)-blastn$jc_dist_round[k]))
}

jc_dist_tibble <- tibble(jc_dist = 0:max(blastn$jc_dist_round), bp = colSums(jc_dist_mat)) %>%
  mutate(jc_dist = jc_dist/100)

ggplot(jc_dist_tibble, aes(jc_dist, bp)) + geom_bar(colour="black", fill="#DD8888", stat="identity")
