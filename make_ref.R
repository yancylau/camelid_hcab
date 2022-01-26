#!/usr/bin/env Rscipt
# 
# Make the germline V gene annotation file (IMGT numbering system)
# 
# Reference: 
# ncbi.github.io/igblast/cook/How-to-set-up.html


library(tidyverse)

ighv <- read_csv("../../refs/ighv.csv") 

ighv <- read_tsv("ighv.tsv", col_names = F) %>% set_names(c("id", "imgt_nt")) 


ndm <- ighv %>% 
  mutate(fwr1 = toupper(str_remove_all(substr(imgt_nt, 1, 78), "\\.")),
         cdr1 = toupper(str_remove_all(substr(imgt_nt, 79, 114), "\\.")),
         fwr2 = toupper(str_remove_all(substr(imgt_nt, 115, 165), "\\.")),
         cdr2 = toupper(str_remove_all(substr(imgt_nt, 166, 195), "\\.")),
         fwr3 = toupper(str_remove_all(substr(imgt_nt, 196, 312), "\\."))) %>%
  mutate(fwr1_len = nchar(fwr1),
         cdr1_len = nchar(cdr1),
         fwr2_len = nchar(fwr2),
         cdr2_len = nchar(cdr2),
         fwr3_len = nchar(fwr3)) %>%
  mutate(fwr1_start = 1, fwr1_end = fwr1_len,
         cdr1_start = fwr1_len + 1, 
         cdr1_end = fwr1_len + cdr1_len,
         fwr2_start = fwr1_len + cdr1_len + 1, 
         fwr2_end = fwr1_len + cdr1_len + fwr2_len,
         cdr2_start = fwr1_len + cdr1_len + fwr2_len + 1, 
         cdr2_end = fwr1_len + cdr1_len + fwr2_len + cdr2_len,
         fwr3_start = fwr1_len + cdr1_len + fwr2_len + cdr2_len + 1, 
         fwr3_end = fwr1_len + cdr1_len + fwr2_len + cdr2_len + fwr3_len) %>%
  select(id, 
         fwr1_start, fwr1_end, 
         cdr1_start, cdr1_end, 
         fwr2_start, fwr2_end, 
         cdr2_start, cdr2_end, 
         fwr3_start, fwr3_end)

write_tsv(ndm, "ndm.imgt", col_names = F)



##
pdm <- ndm %>% 
  mutate(fwr1_aa_end = fwr1_end / 3,
         cdr1_aa_end = cdr1_end / 3,
         fwr2_aa_end = fwr2_end / 3,
         cdr2_aa_end = cdr2_end / 3,
         fwr3_aa_end = fwr3_end / 3,
         fwr1_aa_start = 1,
         cdr1_aa_start = fwr1_aa_end + 1,
         fwr2_aa_start = cdr1_aa_end + 1,
         cdr2_aa_start = fwr2_aa_end + 1,
         fwr3_aa_start = cdr2_aa_end + 1) %>%
  select(id, 
         fwr1_aa_start, fwr1_aa_end, 
         cdr1_aa_start, cdr1_aa_end, 
         fwr2_aa_start, fwr2_aa_end, 
         cdr2_aa_start, cdr2_aa_end, 
         fwr3_aa_start, fwr3_aa_end)

write_tsv(pdm, "pdm.imgt", col_names = F)


# df1 <- ighv %>%
#   mutate(
#     fwr1 = substr(full_sequence, fwr1_start, fwr1_end),
#     cdr1 = substr(full_sequence, cdr1_start, cdr1_end),
#     fwr2 = substr(full_sequence, fwr2_start, fwr2_end),
#     cdr2 = substr(full_sequence, cdr2_start, cdr2_end),
#     fwr3 = substr(full_sequence, fwr3_start, fwr3_end)
#   ) %>%
#   select(sequence_id, fwr1, cdr1, fwr2, cdr2, fwr3)
# 
#   mutate(fwr1_len = nchar(fwr1),
#          cdr1_len = nchar(cdr1),
#          fwr2_len = nchar(fwr2),
#          cdr2_len = nchar(cdr2),
#          fwr3_len = nchar(fwr3)) %>%
#   select(sequence_id, fwr1_len, cdr1_len, fwr2_len, cdr2_len, fwr3_len)

people <- data.frame(city, age, sex)

city <- c('bj', 'sh', 'cd', 'sh', 'bj')
age <- c(12, 12, 14, 15, 40)
sex =c("F", "F", "M", "M", "F")
people <- data.frame(city, age, sex)
people
