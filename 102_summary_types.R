#!/usr/bin/env Rscript
#

load("data/v_gdna/v_gdna.Rda")
load("data/v_gdna/otus.Rda")


seq_types <- otus %>% 
  group_by(sample, type) %>% summarise(count = n()) %>% 
  group_by(sample) %>% mutate(prop = count / sum(count)) %>%
  mutate(lab = paste0(round(prop*100, 2), "% (", count, ")")) 


# Sort samples by VHH proportion
sample_level <- seq_types %>% subset(type == "VHH" & sample != "gDNA-wild-1") %>% arrange(-prop) %>% select(sample) %>% pull()
sample_level <- c( sample_level, "gDNA-wild-1")
seq_types <- seq_types %>% mutate(sample = factor(sample, levels = sample_level))


# Chi square test ---------------------------------------------------------
# By samples: X-squared = 10.459, df = 7, p-value = 0.164
chisq_test_by_sample <- seq_types %>% select(sample, type, count) %>%
  spread(type, count) %>% column_to_rownames("sample") %>% chisq.test()

# By species: X-squared = 1.2938, df = 1, p-value = 0.2554
chisq_test_by_species <- seq_types %>% 
  separate(sample, c(NA, "species", NA), sep = "-") %>% 
  group_by(species, type) %>% summarise(count = sum(count)) %>% 
  spread(type, count) %>% column_to_rownames("species") %>%  chisq.test()


# Plot
p_types_by_samples <- ggplot(seq_types, aes(x = sample, y = count, fill = type)) +
  geom_bar(stat = "identity", position = "fill", width = 0.9) +
  geom_text(aes(y = prop, label = lab), position = position_fill(vjust = 0.5), size = 3, angle = 90) +
  # scale_y_continuous(limits = c(0, 90), breaks = seq(0, 90, 20)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "",
       y = "Count",
       fill = "Type", 
       caption = ) +
  theme(legend.position = "top") 
  
p_types_by_samples

p_type_by_speices <- seq_types %>% select(sample, type, count) %>%
  mutate(species = ifelse(sample =="gDNA-wild-1", "wild", "domestic")) %>%
  group_by(type, species) %>% summarise(count = sum(count)) %>%
  group_by(species) %>% mutate(prop = count/sum(count)) %>%
  mutate(lab = paste0(round(prop*100, 2), "% (", count, ")")) %>%
  ggplot(aes(x = species, y = prop, fill = type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = prop,label = lab), angle = 90, position = position_fill(vjust = 0.5)) +
  theme_classic() 
p_type_by_speices


pdf("results/v_gdna/sequence_types.pdf", width = 6, height = 5)
p_types_by_samples
p_type_by_speices 
dev.off()




