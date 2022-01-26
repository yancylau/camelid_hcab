#!/usr/bin/env Rscript

library(tidytext)
source("R/analyse_cysteins.R")
load("data/camelids/v_imgt_aa.Rda")


v_gdna <- v_imgt_aa %>%
  select(id, family, species, molecule, subgene, type, from, imgt_aa) %>%
  mutate(subgene = ifelse(subgene == "IGHV3", "IGHV3", "non_IGHV3")) %>%
  mutate(imgt_aa = str_sub(imgt_aa, 1, 104)) %>%
  subset(molecule == "gDNA")

v3_gdna <- v %>% subset(subgene == "IGHV3")



# Compare proporation of cysteines --------------------------------------------

# Count cysteines on V sequences
n_cys <- summary_cysteins(ighv3)
write_csv(n_cys, "results/camelids/cystein_summary.csv")

# Chi-square test on proporation of cysteines
chisq_test <- test_cysteins(n_cys)
write_csv(n_cys, "results/camelids/cystein_chi_square_test.csv")

# Plot cystein count
pdf("results/camelids/cystein_proporation.pdf", width = 6, height = 6)
plot_cystein(n_cys)
dev.off()




# Cysteine Positions ------------------------------------------------------
c_pattern <- function(v, sp, tp) {
  v %>% subset(species == sp) %>% subset(type == tp) %>%
    cystein_pattern() %>% 
    mutate(species = sp) %>%
    mutate(type = tp)
}


cys_vhh_cb <- c_pattern(v, "Camelus bactrianus", "VHH")
cys_vhh_cd <- c_pattern(v, "Camelus dromedarius", "VHH")
cys_vhh_vp <- c_pattern(v, "Vicugna pacos", "VHH")
cys_vh_cb <- c_pattern(v, "Camelus bactrianus", "VH")
cys_vh_cd <- c_pattern(v, "Camelus dromedarius", "VH")
cys_vh_vp <- c_pattern(v, "Vicugna pacos", "VH")

cys_patterns <- cys_vhh_cb %>% bind_rows(cys_vhh_cd) %>% bind_rows(cys_vhh_vp) %>%
  bind_rows(cys_vh_cb) %>% bind_rows(cys_vh_cd) %>% bind_rows(cys_vh_vp) %>%
  replace(is.na(.), 0) 

write_csv(cys_patterns, "results/camelids/cys_patterns.csv")



## Cystein distribution --------------------------------------------------------
cys_positions <- summary_c_positions(v)

pdf("results/camelids/cystein_distrubution.pdf", width = 8, height = 6)
plot_cystein_distribution(cys_positions)
dev.off()





# Sub-types of VHH ---------------------------------------------------------
p_pattern <- ggplot(d, aes(x = species, y = prop, fill = species)) +
  geom_bar(stat = "identity", width = 0.8) +
  # geom_text(aes(label = lab)) +
  facet_grid(type + pattern ~ .) +
  coord_flip() +
  theme_classic()
p_pattern

pdf("results/camelids/cystein_pattern.pdf", width = 6, height = 10)
p_pattern
dev.off()
