#!/usr/bin/env Rscript


dir.create("results/camelids", showWarnings = F)
load("data/camelids/v_imgt_aa.Rda")
source("R/analyse_composite.R")


frs <- v_imgt_aa %>% 
  subset(molecule == "gDNA" & family == "camelids") %>%
  mutate(frs = str_c(imgt_fr1, imgt_fr2, imgt_fr3)) %>%
  select(type, species, frs)


# Main function
analysis_diff_composite <- function(sp, path) {
  path <- paste0(path, "/", str_replace_all(tolower(sp), " ", "_"))
  
  # Composite
  frs %>% subset(species == sp) %>%
    composite_plot(paste0(path, "_aa_composite.pdf"), width = 6, height = 8)
  frs %>% subset(species == sp) %>%
    composite_test(paste0(path, "_aa_composite.csv"))
  
  # Chemical
  frs %>% subset(species == sp) %>%
    chemicals_plot(paste0(path, "_aa_chemical.pdf"))
  frs %>% subset(species == sp) %>%
    chemicals_test(paste0(path, "_aa_cehmical.csv"))
}

analysis_diff_composite("Camelus bactrianus", "results/camelids")
analysis_diff_composite("Camelus dromedarius", "results/camelids")
analysis_diff_composite("Vicugna pacos", "results/camelids")


