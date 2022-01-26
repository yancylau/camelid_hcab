#!/usr/bin/env Rscript
#
# summary V mapping to the human germline V
#


load("data/djc_cdna/vdjc.Rda")
library(rstatix)


## D subfamily
d_subfamily <- vdjc %>% 
  group_by(d_call) %>% summarise(count = n()) %>%
  mutate(prop = count / sum(count)) %>%
  mutate(label = paste0(round(prop*100, 2), "% (", count, ")")) %>%
  mutate(gene = "D") %>%  rename(subfamily = d_call)
p_d_subfamily <- ggplot(d_subfamily, aes(x = gene, y = prop, fill = subfamily)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(fill = "V call subfamily")
p_d_subfamily


## J subfamily
j_subfamily <- vdjc %>% 
  group_by(j_call) %>%
  summarise(count = n()) %>%
  mutate(prop = count / sum(count)) %>%
  mutate(label = paste0(round(prop*100, 2), "% (", count, ")")) %>%
  mutate(gene = "J") %>%  rename(subfamily = j_call)
p_j_subfamily <- ggplot(j_subfamily, aes(x = gene, y = prop, fill = subfamily)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(fill = "V call subfamily")
p_j_subfamily




vdjc_identity <- vdjc %>% select(d_identity, j_identity, ch2_identity)  %>%
  gather(key = "gene", value = "identity") 


p_vdjc_identity <- ggplot(vdjc_identity, aes(x = gene, y = identity, color = gene)) +
  geom_violin() +
  # geom_boxplot(width=0.1) +
  # geom_boxplot() +
  # geom_jitter(width = 0.2, size = 0.1) +
  theme_classic() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 0, colour = "black", hjust=0.95,vjust=0.3),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "", y = "Identity", color = "Type")
p_vdjc_identity


pdf("results/djc_cdna/djc_alignment_summary.pdf")
p_d_subfamily
p_j_subfamily
p_vdjc_identity
dev.off()
