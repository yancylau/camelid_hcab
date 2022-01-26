#!/usr/bin/env Rscript
#
# Merge V sequences of Bactrian camels with other species

dir.create("results/camelids", showWarnings = F)
dir.create("data/camelids", showWarnings = F)


load("data/v_gdna/v_gdna.Rda")
load("data/refs/v_imgt_aa_refs.Rda")
# load("data/refs/v.Rda")

v_bactrian <- v_gdna %>% 
  rename(
    imgt_fr1 = fwr1_imgt_aa,
    imgt_fr2 = fwr2_imgt_aa,
    imgt_fr3 = fwr3_imgt_aa,
    imgt_cdr1 = cdr1_imgt_aa,
    imgt_cdr2 = cdr2_imgt_aa
  ) %>%
  mutate(imgt_aa = paste0(imgt_fr1, imgt_cdr1, imgt_fr2, imgt_cdr2, imgt_fr3)) %>%
  select(id, type, imgt_aa, imgt_fr1, imgt_fr2, imgt_fr3, imgt_cdr1, imgt_cdr2) %>%
  mutate(
    subgene = "IGHV3", 
    species = "Camelus bactrianus", 
    family = "camelids",
    molecule = "gDNA",
    from = "this_study"
  )

  
v_imgt_aa <- v_imgt_aa_refs %>% 
  # mutate(sequence_id = paste0("seq", rownames(.))) %>%
  bind_rows(v_bactrian)

# Save metadata to file
save(v_imgt_aa, file = "data/camelids/v_imgt_aa.Rda")




