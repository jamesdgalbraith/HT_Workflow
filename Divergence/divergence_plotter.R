suppressMessages(library(tidyverse))
# species_list <- read_tsv("", col_names = c("species", "genome"))

# Read in RepeatMasker out
rm_out_aipLae <- read_table("~/HT_Workflow/Divergence/assembly_20171114.fasta.out", col_names = c("SW", 'div', "del", "ins", "sseqid", "sstart", "send", "sleft", "strand", "qseqid", "class", "qstart", "qend", "qleft", "ID", "dupl"), skip = 3)
rm_out_emyIji <- read_table("~/HT_Workflow/Divergence/emyIji_1.0.fasta.out", col_names = c("SW", 'div', "del", "ins", "sseqid", "sstart", "send", "sleft", "strand", "qseqid", "class", "qstart", "qend", "qleft", "ID", "dupl"), skip = 3)
rm_out_hydMel <- read_table("~/HT_Workflow/Divergence/hydMel_1.0.fasta.out", col_names = c("SW", 'div', "del", "ins", "sseqid", "sstart", "send", "sleft", "strand", "qseqid", "class", "qstart", "qend", "qleft", "ID", "dupl"), skip = 3)
rm_out_hydCya <- read_table("~/HT_Workflow/Divergence/ASM402372v1.fasta.out", col_names = c("SW", 'div', "del", "ins", "sseqid", "sstart", "send", "sleft", "strand", "qseqid", "class", "qstart", "qend", "qleft", "ID", "dupl"), skip = 3)
rm_out_notScu <- read_table("~/HT_Workflow/Divergence/TS10Xv2-PRI.fasta.out", col_names = c("SW", 'div', "del", "ins", "sseqid", "sstart", "send", "sleft", "strand", "qseqid", "class", "qstart", "qend", "qleft", "ID", "dupl"), skip = 3)
rm_out_pseTex <- read_table("~/HT_Workflow/Divergence/EBS10Xv2-PRI.fasta.out", col_names = c("SW", 'div', "del", "ins", "sseqid", "sstart", "send", "sleft", "strand", "qseqid", "class", "qstart", "qend", "qleft", "ID", "dupl"), skip = 3)

# Filter repeat
rm_out_Proto2_aipLae <- rm_out_aipLae %>%
  filter(qseqid == "Proto2-Snek")
rm_out_Rex1_1_aipLae <- rm_out_aipLae %>%
  filter(qseqid == "Rex1-Snek_1_aipLae")
rm_out_Rex1_2_aipLae <- rm_out_aipLae %>%
  filter(qseqid == "Rex1-Snek_2_aipLae")
rm_out_RTE_aipLae <- rm_out_aipLae %>%
  filter(qseqid == "RTE-Snek_aipLae")

rm_out_Proto2_emyIji <- rm_out_emyIji %>%
  filter(qseqid == "Proto2-Snek")
rm_out_Rex1_1_emyIji <- rm_out_emyIji %>%
  filter(qseqid == "Rex1-Snek_1_aipLae")
rm_out_Rex1_2_emyIji <- rm_out_emyIji %>%
  filter(qseqid == "Rex1-Snek_2_aipLae")
rm_out_RTE_emyIji <- rm_out_emyIji %>%
  filter(qseqid == "RTE-Snek_aipLae")

rm_out_Proto2_hydMel <- rm_out_hydMel %>%
  filter(qseqid == "Proto2-Snek")
rm_out_Rex1_1_hydMel <- rm_out_hydMel %>%
  filter(qseqid == "Rex1-Snek_1_aipLae")
rm_out_Rex1_2_hydMel <- rm_out_hydMel %>%
  filter(qseqid == "Rex1-Snek_2_aipLae")
rm_out_RTE_hydMel <- rm_out_hydMel %>%
  filter(qseqid == "RTE-Snek_aipLae")

rm_out_Proto2_hydCya <- rm_out_hydCya %>%
  filter(qseqid == "Proto2-Snek")
rm_out_Rex1_1_hydCya <- rm_out_hydCya %>%
  filter(qseqid == "Rex1-Snek_1_aipLae")
rm_out_Rex1_2_hydCya <- rm_out_hydCya %>%
  filter(qseqid == "Rex1-Snek_2_aipLae")
rm_out_RTE_hydCya <- rm_out_hydCya %>%
  filter(qseqid == "RTE-Snek_aipLae")

rm_out_Proto2_notScu <- rm_out_notScu %>%
  filter(qseqid == "Proto2-Snek")
rm_out_Rex1_1_notScu <- rm_out_notScu %>%
  filter(qseqid == "Rex1-Snek_1_aipLae")
rm_out_Rex1_2_notScu <- rm_out_notScu %>%
  filter(qseqid == "Rex1-Snek_2_aipLae")
rm_out_RTE_notScu <- rm_out_notScu %>%
  filter(qseqid == "RTE-Snek_aipLae")

rm_out_Proto2_pseTex <- rm_out_pseTex %>%
  filter(qseqid == "Proto2-Snek")
rm_out_Rex1_1_pseTex <- rm_out_pseTex %>%
  filter(qseqid == "Rex1-Snek_1_aipLae")
rm_out_Rex1_2_pseTex <- rm_out_pseTex %>%
  filter(qseqid == "Rex1-Snek_2_aipLae")
rm_out_RTE_pseTex <- rm_out_pseTex %>%
  filter(qseqid == "RTE-Snek_aipLae")

rm_out_Proto2_ophHan <- rm_out_ophHan %>%
  filter(qseqid == "Proto2-Snek")
rm_out_Rex1_1_ophHan <- rm_out_ophHan %>%
  filter(qseqid == "Rex1-Snek_1_aipLae")
rm_out_Rex1_2_ophHan <- rm_out_ophHan %>%
  filter(qseqid == "Rex1-Snek_2_aipLae")
rm_out_RTE_ophHan <- rm_out_ophHan %>%
  filter(qseqid == "RTE-Snek_aipLae")


# plot Rex1_1 repeat
ggplot(rm_out_Rex1_1_aipLae, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 39), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/Rex1_1_aipLae.pdf", width = 10, height = 10, units = "cm")

ggplot(rm_out_Rex1_1_emyIji, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 39), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/Rex1_1_emyIji.pdf", width = 10, height = 10, units = "cm")

ggplot(rm_out_Rex1_1_hydMel, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 39), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/Rex1_1_hydMel.pdf", width = 10, height = 10, units = "cm")

ggplot(rm_out_Rex1_1_hydCya, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 39), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/Rex1_1_hydCya.pdf", width = 10, height = 10, units = "cm")

ggplot(rm_out_Rex1_1_notScu, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 39), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/Rex1_1_notScu.pdf", width = 10, height = 10, units = "cm")

ggplot(rm_out_Rex1_1_pseTex, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 39), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/Rex1_1_pseTex.pdf", width = 10, height = 10, units = "cm")

# plot Rex1_2 repeat
ggplot(rm_out_Rex1_2_aipLae, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 85), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/Rex1_2_aipLae.pdf", width = 10, height = 10, units = "cm")

ggplot(rm_out_Rex1_2_emyIji, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 39), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/Rex1_2_emyIji.pdf", width = 10, height = 10, units = "cm")

ggplot(rm_out_Rex1_2_hydMel, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 39), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/Rex1_2_hydMel.pdf", width = 10, height = 10, units = "cm")

ggplot(rm_out_Rex1_2_hydCya, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 39), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/Rex1_2_hydCya.pdf", width = 10, height = 10, units = "cm")

ggplot(rm_out_Rex1_2_notScu, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 85), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/Rex1_2_notScu.pdf", width = 10, height = 10, units = "cm")

ggplot(rm_out_Rex1_2_pseTex, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 85), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/Rex1_2_pseTex.pdf", width = 10, height = 10, units = "cm")

# plot RTE repeat
ggplot(rm_out_RTE_aipLae, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/RTE_aipLae.pdf", width = 10, height = 10, units = "cm")

ggplot(rm_out_RTE_emyIji, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 39), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/RTE_emyIji.pdf", width = 10, height = 10, units = "cm")

ggplot(rm_out_RTE_hydMel, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 39), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/RTE_hydMel.pdf", width = 10, height = 10, units = "cm")

ggplot(rm_out_RTE_hydCya, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 39), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/RTE_hydCya.pdf", width = 10, height = 10, units = "cm")

ggplot(rm_out_RTE_notScu, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/RTE_notScu.pdf", width = 10, height = 10, units = "cm")

ggplot(rm_out_RTE_pseTex, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/RTE_pseTex.pdf", width = 10, height = 10, units = "cm")

# plot Proto2 repeats
ggplot(rm_out_Proto2_aipLae, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 300), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/Proto2_aipLae.pdf", width = 10, height = 10, units = "cm")

ggplot(rm_out_Proto2_emyIji, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 39), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/Proto2_emyIji.pdf", width = 10, height = 10, units = "cm")

ggplot(rm_out_Proto2_hydMel, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 39), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/Proto2_hydMel.pdf", width = 10, height = 10, units = "cm")

ggplot(rm_out_Proto2_hydCya, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 39), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/Proto2_hydCya.pdf", width = 10, height = 10, units = "cm")

ggplot(rm_out_Proto2_notScu, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 300), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/Proto2_notScu.pdf", width = 10, height = 10, units = "cm")

ggplot(rm_out_Proto2_pseTex, aes(div)) + geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(-1,40), expand = c(0,0), labels = (0:8)*5, breaks = (0:8)*5) +
  scale_y_continuous(limits = c(0, 300), expand = c(0,0)) +
  labs(x = "Percent divergence from consensus", y = "Number of repeats") + theme_bw()

ggsave(filename = "~/HT_Workflow/Divergence/Proto2_pseTex.pdf", width = 10, height = 10, units = "cm")