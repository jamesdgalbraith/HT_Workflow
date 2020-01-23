suppressMessages(library(tidyverse))

# species_list <- read_tsv("", col_names = c("species", "genome"))

# Read in RepeatMasker out
rm_out_aipLae <- read_table("Divergence/assembly_20171114.fasta.out", col_names = c("SW", 'div', "del", "ins", "sseqid", "sstart", "send", "sleft", "strand", "qseqid", "class", "qstart", "qend", "qleft", "ID", "dupl"), skip = 3)
rm_out_notScu <- read_table("Divergence/TS10Xv2-PRI.fasta.out", col_names = c("SW", 'div', "del", "ins", "sseqid", "sstart", "send", "sleft", "strand", "qseqid", "class", "qstart", "qend", "qleft", "ID", "dupl"), skip = 3)
rm_out_pseTex <- read_table("Divergence/EBS10Xv2-PRI.fasta.out", col_names = c("SW", 'div', "del", "ins", "sseqid", "sstart", "send", "sleft", "strand", "qseqid", "class", "qstart", "qend", "qleft", "ID", "dupl"), skip = 3)

rm_out_Rex1_1_aipLae %>%
  filter(div <= 3) %>%
  mutate(length = send - sstart + 1) %>%
  View()


# Filter repeat
rm_out_Proto2_aipLae <- rm_out_aipLae %>%
  filter(qseqid == "Proto2-Snek_aipLae")
rm_out_Proto2_notScu <- rm_out_notScu %>%
  filter(qseqid == "Proto2-Snek_aipLae")
rm_out_Proto2_pseTex <- rm_out_pseTex %>%
  filter(qseqid == "Proto2-Snek_aipLae")

rm_out_Rex1_1_aipLae <- rm_out_aipLae %>%
  filter(qseqid == "Rex1-Snek_1_aipLae")
rm_out_Rex1_1_notScu <- rm_out_notScu %>%
  filter(qseqid == "Rex1-Snek_1_aipLae")
rm_out_Rex1_1_pseTex <- rm_out_pseTex %>%
  filter(qseqid == "Rex1-Snek_1_aipLae")

rm_out_Rex1_2_aipLae <- rm_out_aipLae %>%
  filter(qseqid == "Rex1-Snek_2_aipLae")
rm_out_Rex1_2_notScu <- rm_out_notScu %>%
  filter(qseqid == "Rex1-Snek_2_aipLae")
rm_out_Rex1_2_pseTex <- rm_out_pseTex %>%
  filter(qseqid == "Rex1-Snek_2_aipLae")

rm_out_RTE_aipLae <- rm_out_aipLae %>%
  filter(qseqid == "RTE-Snek_aipLae")
rm_out_RTE_notScu <- rm_out_notScu %>%
  filter(qseqid == "RTE-Snek_aipLae")
rm_out_RTE_pseTex <- rm_out_pseTex %>%
  filter(qseqid == "RTE-Snek_aipLae")

# plot repeat
ggplot(rm_out_Rex1_2_aipLae, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 85), expand = c(0,0))

ggplot(rm_out_Rex1_2_notScu, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 85), expand = c(0,0))

ggplot(rm_out_Rex1_2_pseTex, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 85), expand = c(0,0))
