#!/usr/bin/env Rscript
#
# Plot phylogenetics tree

#  1. Merge v_imgt_nt from bactrian, camelids, and non-camelids
#  2. Remove gaps and convert to fasta
#  3. Run clustal to get aligned fasta
#  4. Get sequence annot from v_imgt_aa files



library(pegas)
library(ape)
library(ggtree)
library(ggnewscale)
library(viridis)


# dir.create("results/v_gdna/phylogenetic", showWarnings = F)



fr4

fa <- character(nrow(fr4) * 2)
fa[c(TRUE, FALSE)] <- paste0(">", fr4$id)
fa[c(FALSE, TRUE)] <- fr4$fr4
writeLines(fa, "data/refs/fr4_aa.fa")


aa <- read.FASTA("data/refs/fr4_aa.fa", type = "AA")


# Calculate distance
phy <- dist.ml(aa) %>% NJ()

# Family
family <- labels %>% select(family)

# Species
species <- fr4 %>% column_to_rownames("id") %>% select(species)

subgene <- fr4 %>% 
  mutate(gene = str_extract(id, "IGHJ[0-9]+")) %>%
  column_to_rownames("id") %>% select(gene)
  

aa <- fr4 %>%
  mutate(a = str_sub(fr4, 6, 6)) %>%
  group_by(species, a) %>%
  mutate(a = ifelse(a %in% c("Q", "R", "K"), "qin", "shu")) %>%
  column_to_rownames("id") %>% select(a) 
  
  


## Visualize tree with ggtree
circ <- ggtree(phy, layout = "rectangular", branch.length = "none", size = 0.2)

gheatmap(
  circ, aa,
  color = NA,
  offset = 0, 
  width = 0.1, 
  colnames = FALSE
) + new_scale_fill() +
  geom_tiplab(size = 1)

p1 <- gheatmap(
  circ, species,
  color = NA,
  offset = 0, 
  width = 0.1, 
  colnames = FALSE
) + new_scale_fill() 
p1


p2 <- gheatmap(
  circ, subgene,
  color = NA,
  offset = 4, 
  width = 0.1, 
  colnames = T
) + 
  # scale_fill_manual(name = "Type") +
  new_scale_fill() +
  geom_tiplab()
p2

pdf("results/djc_cdna/phylogenetic_tree.pdf", height = 20)
p2
dev.off()

p3 <- gheatmap(
  p2, subgene,
  color = NA,
  offset = 8, 
  width = 0.1, 
  colnames = FALSE
) + 
  scale_fill_manual(values = c("#fdae6b", "#fee6ce"), name = "Subgene") +
  # guides(fill = guide_legend(nrow = 3)) +
  new_scale_fill()
# p3


p4 <- gheatmap(
  p3, family,
  color = NA,
  offset = 12, 
  width = 0.1, 
  colnames = FALSE
) + 
  scale_fill_manual(values = c("#0868ac", "#92c5de"), name = "Family") +
  new_scale_fill() +
  guides(fill = guide_legend(ncol = 2)) 
# p4

pdf("results/v_gdna/phylogenetic_tree.pdf")
p4
dev.off()






## Number of cysterin summary --------------------------------------------------

c_groups <- v %>% 
  subset(molecule == "gDNA") %>% select(-molecule) %>%
  rename(full = imgt_aa) %>%
  gather(key = "region", value = "seq", -c(id, family, species, subgene, type)) %>%
  mutate(n_cys = str_count(seq, "C")) %>%
  select(-seq) %>% spread(region, n_cys, fill = 0) %>%
  mutate(pattern = paste0(full, "(", imgt_fr1, "+", imgt_cdr1, "+", imgt_fr2, "+", imgt_cdr2, "+", imgt_fr3, ")")) %>%
  group_by(pattern) %>% add_count(name = "count") %>%
  # mutate(group = ifelse(count > 10, pattern, "others")) %>%
  mutate(group = ifelse(full > 3 || count < 10, "others", pattern)) -> d


# select(id, family, )
# # mutate(n_cysteines = factor(n_cysteines, levels = c("0", "1", "2", "3", "4", "5"))) %>%
# # mutate(region = factor(region, levels = c("full" ,"fwr1", "cdr1", "fwr2", "cdr2", "fwr3"))) %>%
# mutate(label = abs(round(proportion, 2))) %>%
# mutate(label = ifelse(label > 0.1, label, "")) %>%
# mutate(sample = factor(sample, levels = c("B1", "B4", "B5", "B6", "B7", "Y3", "Y4", "Y5", "B4_Y4")))
# 
# select(-id, -imgt_aa) %>% group_by_all() %>% summarise(count = n()) %>%
# group_by(family, species, subgene, type, from) %>%
# mutate(total = sum(count)) %>%
# mutate(prop = count / total) 




# 
# 
# 
# 
# 
# 
# 
# ## Labels
# load("data/v_subgroups.Rda") 
# load("data/merge/vdj.Rda")
# 
# 
# load("data/refs/v_imgt_aa.Rda")
# load("data/v_gdna/v_gdna.Rda")
# 
# 
# type <- vdj %>% subset(molecule=="gDNA") %>% column_to_rownames("id") %>% select(type) 
# species <- vdj %>% subset(molecule=="gDNA") %>% column_to_rownames("id") %>% select(species) 
# subgroup <- v_subgroups %>% 
#   mutate(subgroup = case_when(group == "2(1+0+0+0+1)" ~ "I",
#                               group == "3(1+1+0+0+1)" ~ "II",
#                               group == "3(1+0+1+0+1)" ~ "III",
#                               group == "4(2+1+0+0+1)" ~ "IV",
#                               group == "3(1+0+0+0+2)" ~ "V",
#                               TRUE ~ "Others")) %>%
#   mutate(subgroup = factor(subgroup, levels = c("I", "II", "III", "IV", "V", "Others"))) %>%
#   column_to_rownames("id") %>% select(subgroup) -> d
# 
# ## Labels
# labels <- vdj %>% left_join(v_subgroups, by = c("id", "type", "species")) %>%
#   select(id, type, species, group) %>%
#   column_to_rownames("id")
# 
# #### camelids
# nt <- read.dna("data/fwrs_imgt_nt_gdna.fa", format = "fasta")
# # aa <- trans(nt, 2)
# 
# phy <- dist.dna(
#   nt, 
#   model = "K80", 
#   variance = FALSE,
#   gamma = FALSE, 
#   pairwise.deletion = FALSE, 
#   base.freq = NULL, 
#   as.matrix = FALSE
# ) %>% nj()
# 
# # my_boot <- boot.phylo(
# #   phy, 
# #   nt, 
# #   function(x) nj(dist.gene(x)),
# #   B = 100,
# #   trees = TRUE
# # )
# 
# 
# 
# 
# ## Visualize tree with ggtree
# circ <- ggtree(phy, layout = "circular", branch.length = "none", size = 0.2)
# 
# p1 <- gheatmap(
#   circ, species,
#   color = NA,
#   offset = 0, 
#   width = 0.1, 
#   colnames = FALSE
# ) + 
#   scale_fill_viridis_d(
#     option = "inferno", 
#     name = "Species"
#   ) +
#   new_scale_fill() 
# 
# p2 <- gheatmap(
#   p1, type,
#   color = NA,
#   offset = 4, 
#   width = 0.1, 
#   colnames = FALSE
# ) +
#   scale_fill_manual(values = c("#18908D", "#FFE630"), name = "Type") +
#   new_scale_fill()
# 
# p3 <- gheatmap(
#   p2, subgroup,
#   color = NA,
#   offset = 8, 
#   width = 0.1, 
#   colnames = FALSE
# ) + 
#   # scale_fill_manual(values = c("#d53e4f", "#fc8d59", "#fee08b", "#e6f598", "#99d594", "#3288bd"), name = "Subgroup")
#   scale_fill_viridis_d(option = "D", name = "Subgroup", ) +
#   # theme(legend.position = "bottom") +
#   guides(fill = guide_legend(nrow = 3))
# 
# p3 
# pdf("results/merge/phylogenetic.pdf")
# p3
# dev.off()
# 
# 
# 
# 
