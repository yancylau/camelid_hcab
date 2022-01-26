#!/usr/bin/env Rscript
#


source("R/analyse_usage.R")
load("data/djc_cdna/vdjc.Rda")



# Summary subgene usage
usages_sample <- vdjc %>%
  select(sample, d_call, j_call) %>% 
  gather("gene", "subgene", -sample) %>%
  subset(!is.na(subgene)) %>%
  group_by_all() %>% summarise(count = n()) %>%
  group_by(sample, gene) %>% mutate(prop = count / sum(count)) 


# Gene usage by sample
# Heatmap
usage_d <- usages_sample %>% subset(gene == "d_call")
usage_j <- usages_sample %>% subset(gene == "j_call")
usage_heatmap(usage_d, "results/djc_cdna/heatmap_d_call.pdf")
usage_heatmap(usage_j, "results/djc_cdna/heatmap_j_call.pdf")


# Overall gene usage
# Pie plot
usages_overall <- usages_sample %>% 
  group_by(gene, subgene) %>% summarise(count = sum(count)) %>%
  group_by(gene) %>% mutate(prop = count / sum(count))

pdf("results/djc_cdna/usage_pie_dj.pdf")
usages_overall %>% subset(gene == "d_call") %>% usage_pie()
usages_overall %>% subset(gene == "j_call") %>% usage_pie()
dev.off()


# Average over gene usage 
# Barplot
mean_usages <- usages_sample %>%
  group_by(gene, subgene) %>% summarise(mean = mean(prop), sd = sd(prop)) %>%
  mutate(subgene = fct_reorder(subgene, -mean))

pdf("results/djc_cdna/usage_bar_dj.pdf")
mean_usages %>% subset(gene == "d_call") %>% usage_bar()
mean_usages %>% subset(gene == "d_call") %>% usage_bar()
dev.off()



# IGHJ seqlog logo
ighj <- read_tsv("D:/projects/camel_ig/manuscript/tables/ighj.tsv", col_names = F) %>%
  set_names(c("gene", "sequence"))

pdf("results/djc_cdna/j_seqlogo.pdf", width = 6, height = 2)
ighj %>% subset(gene %in% c("IGHJ4", "IGHJ6")) %>% pull(sequence) %>% seqlogo_plot()
ighj %>% subset(gene %in% c("IGHJ1", "IGHJ2", "IGHJ3", "IGHJ5", "IGHJ7")) %>%
  pull(sequence) %>% seqlogo_plot()
ighj %>% pull(sequence) %>% seqlogo_plot()
dev.off()

