#!/usr/bin/env Rscript
#


dir.create("results/camelids", showWarnings = F)
# load("data/v_gdna/v_gdna.Rda")
# load("data/v_gdna/otus.Rda")
# load("data/refs/v_camelids.Rda")

load("data/camelids/v_imgt_aa.Rda")




# Plot length distribution ----------------------------------------------------
seq_length <- v_imgt_aa %>% 
  subset(family == "camelids" & molecule == "gDNA") %>%
  select(species, type, imgt_fr1, imgt_fr2, imgt_fr3, imgt_cdr1, imgt_cdr2) %>%
  gather(key = "region", value = "seq", -c(species, type)) %>%
  mutate(seq = str_remove_all(seq, "\\.")) %>%
  mutate(length = nchar(seq)) %>%
  group_by(species, type, region, length) %>% summarise(count = n()) %>% 
  arrange(species, region, length, type) %>%
  mutate(region = str_remove_all(region, "imgt_"))

p_cdr2_length <- seq_length %>%
  select(type, length, count) %>%
  mutate(length = as.factor(length)) %>%
  group_by(species, type) %>% mutate(prop = count / sum(count)) %>%
  mutate(label = paste0(round(prop*100, 2), "%\n(", count, ")")) %>%
  mutate(region = str_remove_all(region, "_aa")) %>%
  mutate(region = factor(region, levels = c("fr1", "cdr1", "fr2", "cdr2", "fr3"))) %>%
  mutate(length = factor(length, levels = c(6, 7, 8, 17, 25, 38))) %>%
  subset(region == "cdr2") %>%
  ggplot(aes(x = type, y = prop, fill = length)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  facet_grid(. ~ species) +
  theme_classic() 
p_cdr2_length

pdf("results/camelids/cdr2_length.pdf", width = 8, height = 4)
p_cdr2_length
dev.off()



# Chi square test  --------------------------------------------------------

chisq_test <- v_imgt_aa %>% 
  subset(family == "camelids" & molecule == "gDNA") %>%
  select(species, type, imgt_fr1, imgt_fr2, imgt_fr3, imgt_cdr1, imgt_cdr2) %>%
  gather(key = "region", value = "seq", -c(species, type)) %>%
  mutate(length = nchar(str_remove_all(seq, "\\."))) %>%
  subset(region == "imgt_cdr2") %>%
  select(species, type, length) %>%
  nest_by(species) %>%
  mutate(model = list(chisq.test(data$type, data$length))) %>%
  summarise(broom::tidy(model))


# Results
# species             statistic  p.value parameter method                                                      
# <chr>                   <dbl>    <dbl>     <int> <chr>                                                       
#   1 Camelus bactrianus      54.1  1.89e-13         1 Pearson's Chi-squared test with Yates' continuity correction
# 2 Camelus dromedarius     14.7  6.37e- 4         2 Pearson's Chi-squared test                                  
# 3 Vicugna pacos            3.54 5.98e- 2         1 Pearson's Chi-squared test with Yates' continuity correction


