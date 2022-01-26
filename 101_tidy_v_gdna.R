#!/usr/bin/env Rscript
#
#
# Tidy representative genes (n=170)
#
# 1. Select fields from IgBLAST results
# 2. IMGT numbering for NT and AA sequences
# 3. Determine sequence types based on hallmarks
#          VH      VHH
#  42     V,I      F,Y
#  49       G      E,Q
#  50       L      C,R
#  52     F,W      G,L
# Reference:
# Muyldermans, S. (2013).  Annu. Rev. Biochem. 82, 775–797.

dir.create("data/v_gdna", showWarnings = F)





sequence_id <- read_csv("imported/v_gdna/sequence_id.csv") %>% 
  mutate(otu = as.character(otu))


# Select fiedds -----------------------------------------------------------
igblast <- read_tsv("imported/v_gdna/igblast/atleast3_0.99/igblast.airr")
selected_fields = c("sequence_id", "sequence", "stop_codon", 
                    "v_call", "v_identity",
                    "sequence_alignment", "sequence_alignment_aa", 
                    "fwr1", "cdr1", "fwr2", "cdr2", "fwr3",
                    "fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa", "fwr3_aa")



# IMGT numbering and determine types -----------------------------------------
source("R/imgt_numbering.R")
v_gdna <- igblast %>% select(all_of(selected_fields)) %>%
  mutate(otu = paste0(0:169)) %>% select(-sequence_id) %>% 
  left_join(sequence_id, by = "otu") %>%
  relocate(sequence_id, otu) %>% arrange(sequence_id) %>%
  mutate(
    fwr1_imgt_aa = fr_numbering_aa(fwr1_aa, "FR1"),
    fwr2_imgt_aa = fr_numbering_aa(fwr2_aa, "FR2"),
    fwr3_imgt_aa = fr_numbering_aa(fwr3_aa, "FR3"),
    cdr1_imgt_aa = cdr_numbering_aa(cdr1_aa, "CDR1"),
    cdr2_imgt_aa = cdr_numbering_aa(cdr2_aa, "CDR2")
  ) %>%
  mutate(
    fwr1_imgt = fr_numbering_nt(fwr1, "FR1"),
    fwr2_imgt = fr_numbering_nt(fwr2, "FR2"),
    fwr3_imgt = fr_numbering_nt(fwr3, "FR3"),
    cdr1_imgt = cdr_numbering_nt(cdr1, "CDR1"),
    cdr2_imgt = cdr_numbering_nt(cdr2, "CDR2")
  ) %>%
  unite("imgt_aa", c(fwr1_imgt_aa, cdr1_imgt_aa, fwr2_imgt_aa, cdr2_imgt_aa, fwr3_imgt_aa), remove = FALSE, sep = "") %>%
  unite("imgt_nt", c(fwr1_imgt, cdr1_imgt, fwr2_imgt, cdr2_imgt, fwr3_imgt), remove = FALSE, sep = "") %>% 
  mutate(
    a42 = substr(imgt_aa, 42, 42), 
    a49 = substr(imgt_aa, 49, 49), 
    a50 = substr(imgt_aa, 50, 50), 
    a52 = substr(imgt_aa, 52, 52)
  ) %>%
  mutate(
    type = case_when(
      a49 == "G" ~ "VH",
      a49 %in% c("E", "Q") ~ "VHH",
      a50 == "L" ~ "VH",
      a50 %in% c("C", "R") ~ "VHH",
      TRUE ~ "unknow"
    )
  ) %>%
  rename(id = sequence_id)


write_tsv(v_gdna, "data/v_gdna/v_gdna.tsv")
save(v_gdna, file = "data/v_gdna/v_gdna.Rda")





# http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
#
# IMGT unique numbering (aa/nt)
# FWR1 :     1 -  26 (26)      1-78
# CDR1 :    27 -  38 (12)     81-114
# FWR2 :    39 -  55 (17)    117-165
# CDR2 :    56 -  65 (10)    168-195
# FWR3 :    66 - 104 (39)    198-312
# CDR3 :   105 - 117 (13)    313-351
# FWR4 :   118 - 129（12）   352-287





