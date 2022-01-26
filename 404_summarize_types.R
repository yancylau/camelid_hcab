#!/usr/bin/env Rscript


load("data/v_cdna/v_cdna.Rda")

v_cdna <- v_cdna %>%
  subset(type != "unknown") %>%
  mutate(type = factor(type, levels = c("VH", "VHH")))


# type_summaries <- v_cdna %>% 
#   mutate(type = ifelse(type2 == "traditional", type1, "unknown")) %>%
#   group_by(type, a42, a49, a50, a52) %>% tally() %>%
#   mutate(type = factor(type, levels = c("VH", "VHH", "unknown"))) %>%
#   arrange(type, -n) 
# write_csv(type_summaries, "results/v_cdna/type_summaries.csv")



seq_types <- v_cdna %>% 
  group_by(species, sample, type) %>% summarise(count = n()) %>%
  group_by(species, sample) %>% mutate(prop = count / sum(count)) %>%
  mutate(label = paste0(round(prop, 3) * 100, "%"))


sample_level_wild <- seq_types %>% 
  subset(type == "VHH" & species == "wild") %>% 
  arrange(prop) %>% select(sample) %>% pull()
sample_level_dome <- seq_types %>% 
  subset(type == "VHH" & species == "dome") %>% 
  arrange(prop) %>% select(sample) %>% pull()
sample_level <- c(sample_level_dome, sample_level_wild)


chisq_test <- seq_types %>% ungroup() %>% 
  subset(species == "dome") %>%
  select(sample, type, count)  %>% spread(type, count) %>% 
  column_to_rownames("sample") %>% chisq.test()



p_types_individuals <- seq_types %>%
  mutate(sample = factor(sample, levels = sample_level)) %>%
  ggplot(aes(x = sample, y = prop, fill = type)) +
  geom_bar(stat = "identity", position = position_stack(reverse = T), width = 0.9) +
  # coord_polar(theta = "y")  +
  geom_text(aes(y = prop, label = label), position = position_stack(reverse = T, vjust = 0.5), size = 2.5) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  labs(y = "Proprotion", fill = "Type", caption = paste0("p_value = ", chisq_test$p.value))

# p-value < 2.2e-16

p_types_individuals 


### Summary sequence type by species
chisq_test <- seq_types %>% group_by(species, type) %>% summarise(count = sum(count)) %>%
  spread(type, count) %>%
  column_to_rownames("species") %>%
  chisq.test()


p_type_species <- seq_types %>% 
  group_by(species, type) %>% summarise(count = sum(count)) %>%
  group_by(species) %>% mutate(prop = count/sum(count)) %>%
  mutate(lab = paste0(round(prop, 3) * 100, "%")) %>%
  ggplot(aes(x = species, y = prop, fill = type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label =  lab), position = position_stack(vjust = 0.5, reverse = F)) +
  theme_classic() +
  labs(y = "Proprotion", fill = "Type", caption = paste0("p_value = ", chisq_test$p.value))

# p-value < 2.2e-16

pdf("results/v_cdna/sequence_types.pdf")
p_types_individuals
p_type_species 
dev.off()





