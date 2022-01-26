#!/usr/bin/env Rscript

library(tidytext)
library(pegas)
library(ape)



## Bactrian camels -----------------------------------------------------------
fa <- read.dna("data/v_gdna/v_gdna.fa", format = "fasta")

vh_ids <- labels(fa) %>% data.frame() %>% set_names("id") %>% 
  separate(id, into = c(NA, "type", NA), sep = "_", remove = F) %>%
  subset(type == "vh") %>%
  select(id) %>% pull()

vhh_ids <- labels(fa) %>% data.frame() %>% set_names("id") %>% 
  separate(id, into = c(NA, "type", NA), sep = "_", remove = F) %>%
  subset(type == "vhh") %>%
  select(id) %>% pull()

vh_nt_div <- fa[vh_ids,] %>% nuc.div(variance = FALSE, pairwise.deletion = FALSE)
vhh_nt_div <- fa[vhh_ids,] %>% nuc.div(variance = FALSE, pairwise.deletion = FALSE)


nt_divs_bactrian <- data.frame(
  "species" = "Camelus bactrianus",  
  "subgene" = "IGHV3",
  "type" = c("VH", "VHH"),
  "nt_diversity" = c(vh_nt_div, vhh_nt_div),
  "family" = "camelids",
  "molecule_type" = "gDNA",
  "n_seqs" = c(115, 55)
)





## Camelids ----------------------------------------------
v_fa_camelids <- read.dna("data/refs/v_imgt_nt_camelids.fa", format = "fasta")

id_list_type <- labels(v_fa_camelids) %>% data.frame() %>% set_names("id") %>%
  separate(id, c(NA, "species", "molecule", "type"), sep = "\\|", remove = F) %>%
  group_split(species, molecule, type)

nt_divs_camelids_type <- lapply(
  id_list_type, function(df) {
    n_seqs <- nrow(df)
    sp <- df$species[1]
    molecule_type <- df$molecule[1]
    type <- df$type[1]
    ids <- df %>% select(id) %>% pull()
    
    fa <- v_fa_camelids[ids,]
    #print(sp)
    #print(fa)
    nt_div <- nuc.div(fa, variance = FALSE, pairwise.deletion = FALSE)
    
    data.frame("species" = sp, "n_seqs"= n_seqs, "molecule_type" = molecule_type, 
               "type" = type, "nt_diversity" = nt_div)
  }
) %>% bind_rows() %>% mutate(family = "camelids", subgene = "IGHV3")



id_list <- labels(v_fa_camelids) %>% data.frame() %>% set_names("id") %>%
  separate(id, c(NA, "species", "molecule", NA), sep = "\\|", remove = F) %>%
  group_split(species, molecule)

nt_divs_camelids_all <- lapply(
  id_list, function(df) {
    n_seqs <- nrow(df)
    sp <- df$species[1]
    molecule_type <- df$molecule[1]
    ids <- df %>% select(id) %>% pull()
    
    fa <- v_fa_camelids[ids,]
    #print(sp)
    #print(fa)
    nt_div <- nuc.div(fa, variance = FALSE, pairwise.deletion = FALSE)
    
    data.frame("species" = sp, "n_seqs"= n_seqs, "molecule_type" = molecule_type, "nt_diversity" = nt_div)
  }
) %>% bind_rows() %>% mutate(family = "camelids", subgene = "IGHV3", type = "all")


nt_divs_camelids <- bind_rows(nt_divs_camelids_all, nt_divs_camelids_type)





## Nucleotide diveristy  of non-camelids ----------------------------------------
fasta_seqs <- read.dna("data/refs/v_imgt_nt_non_camelids_truncated.fa", format = "fasta")

id_list <- labels(fasta_seqs) %>% data.frame() %>% set_names("id") %>%
  separate(id, into = c(NA, "subgene", "species", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA), sep = "\\|", remove = F) %>%
  mutate(subgene = str_extract(subgene, "IGHV[0-9]+")) %>%
  separate(species, into = "species", sep = "_", extra = "drop") %>%
  mutate(subgene = ifelse(subgene == "IGHV3", "IGHV3", "non_IGHV3")) %>%
  group_split(species, subgene)

# IGHV3 non_IGHV3
# Bos taurus                   0        15
# Canis lupus familiaris      35         2
# Danio rerio                  1        39
# Equus caballus               3        31
# Gallus gallus                0         2
# Homo sapiens               136       147
# Macaca fascicularis         36        30
# Macaca mulatta              69        56
# Mus musculus                16       242
# Oncorhynchus mykiss          2        58
# Ornithorhynchus anatinus     4        32
# Oryctolagus cuniculus        0        38
# Rattus norvegicus            5       146
# Salmo salar                  3        88
# Sus scrofa                   0        15
# Vicugna pacos                5         2


nt_divs_non_camelids <- lapply(
  id_list, function(df) {
    n_seqs <- nrow(df)
    sp <- df$species[1]
    subgene <- df$subgene[1]
    ids <- df %>% select(id) %>% pull()
    
    fa <- fasta_seqs[ids]
    
    #print(sp)
    #print(fa)
    nt_div <- nuc.div(fa, variance = FALSE, pairwise.deletion = FALSE)
    
    data.frame("species" = sp, "subgene" = subgene, "n_seqs"= n_seqs, "nt_diversity" = nt_div)
  }
) %>% bind_rows() %>% mutate(family = "non-camelids",  type = "VH", molecule_type = "gDNA")




nt_divs <- bind_rows(nt_divs_bactrian, nt_divs_camelids) %>%
  bind_rows(nt_divs_non_camelids) %>%
  subset(!is.na(nt_diversity)) %>%
  mutate(lab = paste0(round(nt_diversity, 3), " (n=", n_seqs, ")")) 

write_csv(nt_divs, "results/nelceotide_diversity.csv")



## Plot diversity
p_nt_diversity_all <- nt_divs %>%
  # subset(n_seqs >= 5) %>%
  subset(molecule_type == "gDNA") %>% subset(type != "all") %>%
  subset(species != "Vicugna pacos") %>%
  mutate(type = factor(type, levels = c("VHH", "VH"))) %>%
  # ggplot(aes(x = reorder(species, -nt_diversity), y = nt_diversity)) +
  # scale_x_reordered() +
  ggplot(aes(x = reorder_within(species, nt_diversity, list(type, subgene)), y = nt_diversity)) +
  scale_x_reordered() +
  geom_bar(stat = "identity", aes(fill = species), width = 0.8) +
  geom_text(aes(y = nt_diversity + 0.05, label = lab)) +
  facet_grid(family + subgene + type ~ ., space = "free", scales = "free") +
  theme_classic() + 
  # theme(axis.text.x = element_text(angle = 90, colour = "black", hjust = 0.95, vjust = 0.2)) +
  labs(x = "", y = "Nucleotide diversity", fill = "Species") +
  guides(fill = guide_legend(ncol = 1)) +
  coord_flip() 

p_nt_diversity_all



p_nt_diversity <- nt_divs %>%
  subset(n_seqs >= 5) %>%
  subset(molecule_type == "gDNA") %>% subset(type != "all") %>%
  subset(species != "Vicugna pacos") %>%
  mutate(type = factor(type, levels = c("VHH", "VH"))) %>%
  # ggplot(aes(x = reorder(species, -nt_diversity), y = nt_diversity)) +
  # scale_x_reordered() +
  ggplot(aes(x = reorder_within(species, nt_diversity, list(type, subgene)), y = nt_diversity)) +
  scale_x_reordered() +
  geom_bar(stat = "identity", aes(fill = species), width = 0.8) +
  geom_text(aes(y = nt_diversity + 0.05, label = lab), size = 3) +
  facet_grid(family + subgene + type ~ ., space = "free", scales = "free") +
  theme_classic() + 
  # theme(axis.text.x = element_text(angle = 90, colour = "black", hjust = 0.95, vjust = 0.2)) +
  labs(x = "", y = "Nucleotide diversity", fill = "Species") +
  guides(fill = guide_legend(ncol = 1)) +
  coord_flip() 

p_nt_diversity


pdf("results/v_gdna/nucleotide_diversity.pdf",  width = 10, height = 8)
p_nt_diversity
p_nt_diversity_all
dev.off()

