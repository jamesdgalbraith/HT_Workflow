library(tidyverse)
library(plyranges)
library(BSgenome)
library(gridExtra)
library(grid)
library(svglite)
#perform the blast
blast=read.table(text=system("blastn -query ~/Analysis/Snake/HT_Workflow/Proto2-Snek/Proto2-Snek.fasta -db ~/Analysis/Genomes/Aipysurus_laevis/assembly_20171114.fasta -evalue 1e-07 -outfmt 6", intern = TRUE)) %>%
  as_tibble() %>%
  filter(V3 >= 97) %>%
  mutate(V3 = 100-V3)

#TE consensus size
cons_len=system("grep -v '>'  ~/Analysis/Snake/HT_Workflow/Proto2-Snek/Proto2-Snek.fasta  | wc | awk '{print $3-$1}'", intern = TRUE)
#list of almost full length fragments
full=blast %>% filter(V4 >= as.double(cons_len) * 0.9)

# create coverage tibble
coverage=matrix(rep(0, length(blast$V1)*as.numeric(cons_len)), byrow = T, ncol = as.numeric(cons_len))
for(i in 1:length(blast$V1)){
  coverage[i,]<-c(rep(0,blast$V7[i]-1),rep(1,abs(blast$V8[i]-blast$V7[i])+1), rep(0,as.numeric(cons_len)-blast$V8[i]))
}
coverage <- tibble(y = colSums(coverage))
coverage$x <- 0:(nrow(coverage)-1)

p1 <- ggplot() + geom_segment(aes(x=blast$V8, xend=blast$V7, y = blast$V3, yend = blast$V3), colour = "grey43", size = 0.5) +
  geom_segment(aes(x=full$V8, xend=full$V7, y = full$V3, yend = full$V3), colour = "red", size = 1) + geom_line(mapping = aes(x = coverage$x, y = coverage$y/max(coverage$y) * 3), colour = "blue") +
  scale_y_continuous(expand = c(0,0), name = "Divergence from consensus", lim = c(0,3.1), sec.axis = sec_axis(~(./max(blast$V3)*max(coverage$y)), name = "Consensus coverage (bp)")) +
  scale_x_continuous(name = "Consensus sequence", expand = c(0,0)) + theme_bw() +
  labs(title = "Proto2-Snek") + theme(plot.title = element_text(hjust = 0.5))
p1
ggsave(filename = "coverage_plots/Proto2-Snek_coverage.svg", plot = p1, width = 100, height = 100, units = "mm")


##### RTE-Snek
#perform the blast
blast=read.table(text=system("blastn -query ~/Analysis/Snake/HT_Workflow/RTE-Snek/RTE-Snek.fasta -db ~/Analysis/Genomes/Aipysurus_laevis/assembly_20171114.fasta -evalue 1e-07 -outfmt 6", intern = TRUE)) %>%
  as_tibble() %>%
  filter(V3 >= 97) %>%
  mutate(V3 = 100-V3)

#TE consensus size
cons_len=system("grep -v '>'  ~/Analysis/Snake/HT_Workflow/RTE-Snek/RTE-Snek.fasta  | wc | awk '{print $3-$1}'", intern = TRUE)
#list of almost full length fragments
full=blast %>% filter(V4 >= as.double(cons_len) * 0.9)

# create coverage tibble
coverage=matrix(rep(0, length(blast$V1)*as.numeric(cons_len)), byrow = T, ncol = as.numeric(cons_len))
for(i in 1:length(blast$V1)){
  coverage[i,]<-c(rep(0,blast$V7[i]-1),rep(1,abs(blast$V8[i]-blast$V7[i])+1), rep(0,as.numeric(cons_len)-blast$V8[i]))
}
coverage <- tibble(y = colSums(coverage))
coverage$x <- 0:(nrow(coverage)-1)

p2 <- ggplot() + geom_segment(aes(x=blast$V8, xend=blast$V7, y = blast$V3, yend = blast$V3), colour = "grey43", size = 0.5) +
  geom_segment(aes(x=full$V8, xend=full$V7, y = full$V3, yend = full$V3), colour = "red", size = 1) + geom_line(mapping = aes(x = coverage$x, y = coverage$y/max(coverage$y) * 3), colour = "blue") +
  scale_y_continuous(expand = c(0,0), name = "Divergence from consensus", lim = c(0,3.1), sec.axis = sec_axis(~(./max(blast$V3)*max(coverage$y)), name = "Consensus coverage (bp)")) +
  scale_x_continuous(name = "Consensus sequence", expand = c(0,0)) + theme_bw() +
  labs(title = "RTE-Snek") + theme(plot.title = element_text(hjust = 0.5))
p2
ggsave(filename = "coverage_plots/RTE-Snek_coverage.svg", plot = p2, width = 100, height = 100, units = "mm")

##### Rex1-Snek_1
#perform the blast
blast=read.table(text=system("blastn -query ~/Analysis/Snake/HT_Workflow/Rex1-Snek_1/Rex1-Snek_1.fasta -db ~/Analysis/Genomes/Aipysurus_laevis/assembly_20171114.fasta -evalue 1e-07 -outfmt 6", intern = TRUE)) %>%
  as_tibble() %>%
  filter(V3 >= 97) %>%
  mutate(V3 = 100-V3)

#TE consensus size
cons_len=system("grep -v '>'  ~/Analysis/Snake/HT_Workflow/Rex1-Snek_1/Rex1-Snek_1.fasta  | wc | awk '{print $3-$1}'", intern = TRUE)
#list of almost full length fragments
full=blast %>% filter(V4 >= as.double(cons_len) * 0.9)

# create coverage tibble
coverage=matrix(rep(0, length(blast$V1)*as.numeric(cons_len)), byrow = T, ncol = as.numeric(cons_len))
for(i in 1:length(blast$V1)){
  coverage[i,]<-c(rep(0,blast$V7[i]-1),rep(1,abs(blast$V8[i]-blast$V7[i])+1), rep(0,as.numeric(cons_len)-blast$V8[i]))
}
coverage <- tibble(y = colSums(coverage))
coverage$x <- 0:(nrow(coverage)-1)

p3 <- ggplot() + geom_segment(aes(x=blast$V8, xend=blast$V7, y = blast$V3, yend = blast$V3), colour = "grey43", size = 0.5) +
  geom_segment(aes(x=full$V8, xend=full$V7, y = full$V3, yend = full$V3), colour = "red", size = 1) + geom_line(mapping = aes(x = coverage$x, y = coverage$y/max(coverage$y) * 3), colour = "blue") +
  scale_y_continuous(expand = c(0,0), name = "Divergence from consensus", lim = c(0,3.1), sec.axis = sec_axis(~(./max(blast$V3)*max(coverage$y)), name = "Consensus coverage (bp)")) +
  scale_x_continuous(name = "Consensus sequence", expand = c(0,0)) + theme_bw() +
  labs(title = "Rex1-Snek_1") + theme(plot.title = element_text(hjust = 0.5))
p3
ggsave(filename = "coverage_plots/Rex1-Snek_1_coverage.svg", plot = p3, width = 100, height = 100, units = "mm")

##### Rex1-Snek_2
#perform the blast
blast=read.table(text=system("blastn -query ~/Analysis/Snake/HT_Workflow/Rex1-Snek_2/Rex1-Snek_2.fasta -db ~/Analysis/Genomes/Aipysurus_laevis/assembly_20171114.fasta -evalue 1e-07 -outfmt 6", intern = TRUE)) %>%
  as_tibble() %>%
  filter(V3 >= 97) %>%
  mutate(V3 = 100-V3)

#TE consensus size
cons_len=system("grep -v '>'  ~/Analysis/Snake/HT_Workflow/Rex1-Snek_2/Rex1-Snek_2.fasta  | wc | awk '{print $3-$1}'", intern = TRUE)
#list of almost full length fragments
full=blast %>% filter(V4 >= as.double(cons_len) * 0.9)

# create coverage tibble
coverage=matrix(rep(0, length(blast$V1)*as.numeric(cons_len)), byrow = T, ncol = as.numeric(cons_len))
for(i in 1:length(blast$V1)){
  coverage[i,]<-c(rep(0,blast$V7[i]-1),rep(1,abs(blast$V8[i]-blast$V7[i])+1), rep(0,as.numeric(cons_len)-blast$V8[i]))
}
coverage <- tibble(y = colSums(coverage))
coverage$x <- 0:(nrow(coverage)-1)

p4 <- ggplot() + geom_segment(aes(x=blast$V8, xend=blast$V7, y = blast$V3, yend = blast$V3), colour = "grey43", size = 0.5) +
  geom_segment(aes(x=full$V8, xend=full$V7, y = full$V3, yend = full$V3), colour = "red", size = 1) + geom_line(mapping = aes(x = coverage$x, y = coverage$y/max(coverage$y) * 3), colour = "blue") +
  scale_y_continuous(expand = c(0,0), name = "Divergence from consensus", lim = c(0,3.1), sec.axis = sec_axis(~(./max(blast$V3)*max(coverage$y)), name = "Consensus coverage (bp)")) +
  scale_x_continuous(name = "Consensus sequence", expand = c(0,0)) + theme_bw() +
  labs(title = "Rex1-Snek_2") + theme(plot.title = element_text(hjust = 0.5))

p4
ggsave(filename = "coverage_plots/Rex1-Snek_2_coverage.svg", plot = p4, width = 100, height = 100, units = "mm")

