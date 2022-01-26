

# Convert gDNA sequences to FASTA
imgt_nt_gdna <- vdj %>% subset(molecule == "gDNA") %>% select(id, imgt_nt) %>%
  mutate(imgt_nt = str_replace_all(imgt_nt, "[.]", "-"))
imgt_fa <- character(nrow(imgt_nt_gdna) * 2)
imgt_fa[c(TRUE, FALSE)] <- paste0(">", imgt_nt_gdna$id)
imgt_fa[c(FALSE, TRUE)] <- imgt_nt_gdna$imgt_nt
writeLines(imgt_fa, "data/imgt_nt_gdna.fa")

# imgt_gDNA - FWRs only
imgt_nt_fwrs <- vdj %>% subset(molecule == "gDNA") %>% 
  mutate(imgt_nt_fwr1 = substr(imgt_nt, 1, 78),
         imgt_nt_cdr1 = substr(imgt_nt, 79, 114),
         imgt_nt_fwr2 = substr(imgt_nt, 115, 165),
         imgt_nt_cdr2 = substr(imgt_nt, 166, 195),
         imgt_nt_fwr3 = substr(imgt_nt, 196, 312)) %>%
  mutate(imgt_nt_fwrs = paste0(imgt_nt_fwr1, imgt_nt_fwr2, imgt_nt_fwr3)) %>%
  select(id, imgt_nt_fwrs) %>%
  mutate(imgt_nt_fwrs = str_replace_all(imgt_nt_fwrs, "[.]", "-"))
imgt_fa <- character(nrow(imgt_nt_fwrs) * 2)
imgt_fa[c(TRUE, FALSE)] <- paste0(">", imgt_nt_fwrs$id)
imgt_fa[c(FALSE, TRUE)] <- imgt_nt_fwrs$imgt_nt_fwrs
writeLines(imgt_fa, "data/fwrs_imgt_nt_gdna.fa")


# nt-gDNA
nt_gdna <- vdj %>% subset(molecule == "gDNA") %>% select(id, nt)
fa <- character(nrow(nt_gdna) * 2)
fa[c(TRUE, FALSE)] <- paste0(">", nt_gdna$id)
fa[c(FALSE, TRUE)] <- nt_gdna$nt
writeLines(fa, "data/nt_gdna.fa")


# NT, FRs only
fwrs_nt <- vdj %>% subset(molecule == "gDNA") %>% 
  mutate(nt_fwr1 = substr(imgt_nt, 1, 78),
         nt_cdr1 = substr(imgt_nt, 79, 114),
         nt_fwr2 = substr(imgt_nt, 115, 165),
         nt_cdr2 = substr(imgt_nt, 166, 195),
         nt_fwr3 = substr(imgt_nt, 196, 312)) %>%
  mutate(imgt_fwrs = paste0(nt_fwr1, nt_fwr2, nt_fwr3)) %>%
  mutate(nt_fwrs = str_replace_all(imgt_fwrs, "[-.]", "")) %>%
  select(id, nt_fwrs)
fa <- character(nrow(fwrs_nt) * 2)
fa[c(TRUE, FALSE)] <- paste0(">", nt_gdna$id)
fa[c(FALSE, TRUE)] <- fwrs_nt$nt_fwrs
writeLines(fa, "data/frs_nt_gdna.fa")
