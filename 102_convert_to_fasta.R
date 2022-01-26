


# Convert table to fasta ------------------------------------------------------

load("data/v_gdna/v_gdna.Rda")

convert2fa <- function(imgt_nt) {
  fa <- character(nrow(imgt_nt) * 2)
  fa[c(TRUE, FALSE)] <- paste0(">", imgt_nt$id)
  fa[c(FALSE, TRUE)] <- imgt_nt$imgt_nt
  
  return(fa)
}

v_gdna %>% 
  mutate(id = sequence_id) %>% 
  mutate(imgt_nt = str_replace_all(imgt_nt, "\\.", "-")) %>%
  convert2fa() %>%
  writeLines("data/v_gdna/v_gdna.fa")
