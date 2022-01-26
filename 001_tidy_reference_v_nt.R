#!/usr/bin/env Rscript
#
# Tidy reference sequences
#   1. Camelids from previous studies (GENBank)
#   2. Non-camelids from IMGT-GENE-DB (remove Vicugna pacos): 
#      functional, no partial


dir.create("data/refs", showWarnings = F)
dir.create("results/refs", showWarnings = F)




## IMGT VH/VHH sequence of camelids from previous studies ------------------------

v_imgt_nt_camelids_raw <- read_tsv("imported/camelids/imgt_nt.tsv", col_names = F) %>% select(1:2) %>% set_names("id", "imgt_nt")
id_list <- read_csv("../refs/ighv_camelids/imgt_v_camelids_list.csv")

v_imgt_nt_camelids <- id_list %>% subset(drop == F) %>%
  mutate(type = type_adjusted) %>%
  select(id, type) %>%
  left_join(v_imgt_nt_camelids_raw ) %>%
  unite("id", c(id, type), sep = "|")

save(v_imgt_nt_camelids, file = "data/refs/v_imgt_nt_camelids.Rda")

# > table(imgt_nt_camelids$species, imgt_nt_camelids$molecule_type)
# 
# cDNA gDNA
# Camelus_dromedarius    0   92
# Lama_glama           155    0
# Vicugna_pacos          0   60


## Convert to fasta
convert2fa <- function(imgt_nt) {
  fa <- character(nrow(imgt_nt) * 2)
  fa[c(TRUE, FALSE)] <- paste0(">", imgt_nt$id)
  fa[c(FALSE, TRUE)] <- imgt_nt$imgt_nt
  
  return(fa)
}

convert2fa(v_imgt_nt_camelids) %>% writeLines("data/refs/v_imgt_nt_camelids.fa")




## Non-camelids VH sequence from IMGT GENE DB --------------------------------

# Filter criteron:
#  1. functional: F, [F], (F)
#  2. partial in 5/3 end

imgt_nt_gene_db <- read_tsv("../refs/IMGT_GENE_DB/IGHV/imgt_nt.tsv", col_names = F) %>% set_names("id", "imgt_nt")

v_imgt_nt_non_camelids <- imgt_nt_gene_db %>%
  separate(id, into = c(NA, "subgene", "species", "functionality", NA, NA, NA, NA, NA, NA, NA, NA, "length", "partial", NA, NA), sep = "\\|", remove = F) %>%
  mutate(subgene = str_extract(subgene, "IGHV[0-9]+")) %>%
  separate(species, into = "species", sep = "_", extra = "drop")  %>%
  separate(length, into = c("length", "gaps", "full_length"), sep = "[+=]") %>%
  subset(functionality %in% c("F", "(F)", "[F]")) %>%
  subset(partial == " ") %>%
  mutate(imgt_nt = str_sub(imgt_nt, 1, full_length)) %>%
  mutate(subgene = ifelse(subgene == "IGHV3", "IGHV3", "non_IGHV3")) %>%
  group_by(species) %>% add_count(name = "total") %>% subset(total > 5) %>%
  subset(species != "Vicugna pacos") %>%
  select(species, id, imgt_nt)

save(v_imgt_nt_non_camelids, file = "data/refs/v_imgt_nt_non_camelids.Rda")
write_csv(v_imgt_nt_non_camelids, file = "data/refs/v_imgt_nt_non_camelids.csv")


# Convert table to fasta
# convert2fa(v_imgt_nt_non_camelids) %>% writeLines("data/refs/v_imgt_nt_gene_db.fa")




## Truncate sequences to FR3 end (104-C)
v_imgt_nt_non_camelids_truncated <- read_csv("data/refs/v_imgt_nt_non_camelids_truncated.csv") 

convert2fa <- function(imgt_nt) {
  fa <- character(nrow(imgt_nt) * 2)
  fa[c(TRUE, FALSE)] <- paste0(">", imgt_nt$id)
  fa[c(FALSE, TRUE)] <- imgt_nt$imgt_nt
  
  return(fa)
}
convert2fa(v_imgt_nt_non_camelids_truncated) %>% writeLines("data/refs/v_imgt_nt_non_camelids_truncated.fa")







