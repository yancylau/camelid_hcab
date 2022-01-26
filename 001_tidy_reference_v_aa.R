#!/usr/bin/env Rscript
#
# Tidy references:
#  Sequences of camelids were collected from previous studies (GenBank and papers)
#  Sequences of non-camelids were downloaded from IMGT-GENE/DB


dir.create("data/refs", showWarnings = F)
dir.create("results/refs", showWarnings = F)



# IMGT VH/VHH sequences of camelids -------------------------------------------
v_reference_extend <- read_csv("imported/refs/v_reference_extend.csv")

v_imgt_aa_camelids <- v_reference_extend %>% rename_all(tolower) %>%
  # mutate(type = type_adjusted) %>%
  select(id, from, species, molecule,type, full, fr1, cdr1, fr2, cdr2, fr3) %>%
  rename(imgt_fr1 = fr1, imgt_fr2 = fr2, imgt_fr3 = fr3, imgt_cdr1 = cdr1, imgt_cdr2 = cdr2, imgt_aa = full) %>%
  mutate(subgene = "IGHV3", family = "camelids")
  
save(v_imgt_aa_camelids, file = "data/refs/v_imgt_aa_camelids.Rda")
write_csv(v_imgt_aa_camelids , "results/refs/v_imgt_aa_camelids.csv")



# IMGT VH sequences of non-camelids --------------------------------------------
imgt_aa <- read_tsv("refs/IMGT_GENE_DB/IGHV/imgt_aa.tsv", col_names = F) %>%
  select(1:2) %>% set_names("id", "imgt_aa")

# inserts <- read_csv("../refs/IMGT_GENE_DB/v_imgt_aa_non_camelids_inserts_positions.csv")

v_imgt_aa_non_camelids_raw <- imgt_aa %>%
  separate(id, into = c(NA, "subgene", "species", "functionality", NA, NA, NA, NA, NA, NA, NA, NA, "length", "partial", NA, NA), sep = "\\|", remove = F) %>%
  mutate(subgene = str_extract(subgene, "IGHV[0-9]+")) %>%
  separate(species, into = "species", sep = "_", extra = "drop") %>%
  separate(length, into = c("length", "gaps", "full_length"), sep = "[+=]") %>%
  subset(functionality %in% c("F", "(F)", "[F]")) %>%
  subset(partial == " ") %>%
  select(id, subgene, species, imgt_aa, full_length)

write_csv(v_imgt_aa_non_camelids_raw,  "data/refs/v_imgt_aa_non_camelids_raw.csv")


## Manually remove inserts
v_imft_aa_non_camelids_tidy <- read_csv("refs/IMGT_GENE_DB/v_imft_aa_non_camelids_tidy.csv") 

v_imgt_aa_non_camelids <- v_imft_aa_non_camelids_tidy %>%
  select(id, subgene, species, imgt_aa) %>%
    mutate(
      imgt_fr1 = str_sub(imgt_aa, 1, 26),
      imgt_fr2 = str_sub(imgt_aa, 39, 55),
      imgt_fr3 = str_sub(imgt_aa, 66, 104),
      imgt_cdr1 = str_sub(imgt_aa, 27, 38),
      imgt_cdr2 = str_sub(imgt_aa, 56, 65)
    ) %>%
  mutate(molecule = "gDNA", type = "VH", family = "non_camelids", from = "IMGT_GENE_DB") %>%
  group_by(species) %>% add_count(name = "total") %>% subset(total > 5) %>% select(-total) %>%  
  subset(species != "Vicugna pacos")

save(v_imgt_aa_non_camelids, file = "data/refs/v_imgt_aa_non_camelids.Rda")
write_csv(v_imgt_aa_non_camelids, "data/refs/v_imgt_aa_non_camelids.csv")




# Merge VH sequences of camelids with non-camelids ----------------------------
v_imgt_aa_refs <- bind_rows(v_imgt_aa_camelids, v_imgt_aa_non_camelids)

save(v_imgt_aa_refs, file ="data/refs/v_imgt_aa_refs.Rda")
write_csv(v_imgt_aa_refs, "results/refs/v_imgt_aa_refs.csv")

v_imgt_aa_summary <- v_imgt_aa_refs %>% 
  mutate(subgene = ifelse(subgene == "IGHV3", subgene, "non-IGHV3")) %>%
  group_by(family, subgene, species, type, molecule, from) %>% summarise(count = n()) %>%
  # spread(type, count, fill = 0) %>% mutate(total = VH + VHH) %>%
  select(family, species, subgene, molecule, type, count, from) %>%
  arrange(family, species, subgene)

write_csv(v_summary, "results/refs/v_imgt_aa_refs_summary.csv")






# ## Hallmarks -------------------------------------------------------------------
# 
# # Non-camelids: V3 and non-V3
# hallmarks_non_camelids <- v_imgt_aa_non_camelids %>% 
#   mutate(a42 = substr(imgt_aa, 42, 42), 
#          a49 = substr(imgt_aa, 49, 49), 
#          a50 = substr(imgt_aa, 50, 50), 
#          a52 = substr(imgt_aa, 52, 52)) %>% 
#   mutate(motif = str_c(a42, a49, a50, a52)) %>%
#   mutate(type = "VH", family = "non-camelids") %>%
#   select(family, subgene, type, motif, a42, a49, a50, a52)
# 
# # Camelids: VH and VHH
# hallmarks_camelids <- v_imgt_aa_camelids %>% 
#   rename_all(tolower) %>%
#   mutate(a42 = substr(imgt_aa, 42, 42), 
#          a49 = substr(imgt_aa, 49, 49), 
#          a50 = substr(imgt_aa, 50, 50), 
#          a52 = substr(imgt_aa, 52, 52)) %>% 
#   mutate(motif = str_c(a42, a49, a50, a52)) %>%
#   mutate(subgene = "IGHV3") %>%
#   mutate(family = "camelids") %>%
#   select(family, subgene, type, motif, a42, a49, a50, a52)
# 
# hallmarks <- bind_rows(hallmarks_camelids, hallmarks_non_camelids)
# 
# 
# save(hallmarks, file = "data/refs/hallmarks.Rda")
# write_csv(hallmarks, "results/refs/hallmarks.csv")



