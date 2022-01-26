
load("data/djc_cdna/vdjc.Rda")

v_reference_extend <- read_csv("imported/refs/v_reference_extend.csv")

cdr3_camelids <-  v_reference_extend %>% rename_all(tolower) %>%
  mutate(type = type_adjusted) %>%
  select(id, from, species, molecule,type, cdr3) %>%
  mutate(subgene = "IGHV3", family = "camelids") %>% 
  subset(molecule=="cDNA") %>%
  select(species, type, cdr3) %>%
  subset(species != "Camelus dromedarius") %>%
  subset(!is.na(cdr3)) %>%
  mutate(seq = str_remove_all(cdr3, "\\.")) %>%
  mutate(cdr3_length = nchar(seq)) %>%
  select(species, type, cdr3_length)

cdr3_bactrian <- vdjc %>% select(species,  cdr3_length) %>%
  mutate(type = "VHH", species = "Camelus bactrianus") 

  
p_length_density <- bind_rows(cdr3_bactrian, cdr3_camelids) %>%
  ggplot(aes(x = cdr3_length, group = species)) +
  # geom_bar(stat = "identity", aes(fill = sample)) +
  geom_density(aes(y = ..density.., fill = species), lwd = 0.5, alpha = 0.5, adjust = 2) +
  facet_grid(type ~ .) +
  theme_classic() +
  labs(fill = "Sample", color = "Sample", x = "Length", y = "Porportion")
p_length_density 


pdf("results/djc_cdna/cdr3_length_3.pdf")
p_length_density 
dev.off()
