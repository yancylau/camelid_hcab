#!/usr/bin/env Rscript
#
# summary V mapping to the human germline V
#

library(rstatix)

load("data/v_cdna/v_cdna.Rda")
igblast_human <- read_tsv("imported/v_cdna/igblast/human/atleast3.airr") 



seq_types <- v_cdna %>% select(sequence_id, type)

v_alignment <- igblast_human %>% 
  mutate(sequence_id = str_remove_all(sequence_id, "SAMPLE=|DUPCOUNT=")) %>% 
  separate(sequence_id, c("sequence_id", "sample", "dupcount"), sep = ";") %>%
  right_join(seq_types, by = "sequence_id") %>%
  select(type, v_call, v_identity) %>% 
  separate(v_call, into = "v_call", sep = "\\*", extra = "drop")



v_subfamily <- v_alignment %>% 
  group_by(type, v_call) %>% summarise(count = n()) %>%
  mutate(prop = count / sum(count)) %>%
  mutate(label = paste0(round(prop*100, 2), "% (", count, ")"))



p_v_subfamily <- ggplot(v_subfamily, aes(x = "", y = prop, fill = v_call)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  facet_grid(. ~ type) +
  coord_polar(theta = "y") +
  theme_void() +
  labs(fill = "V call subfamily")
p_v_subfamily



v_identity <- v_alignment %>% select(type, v_identity) 
stat_test <- v_identity %>% t_test(v_identity ~ type) %>% add_xy_position(x = "type")
v_identity %>% subset(type == "VHH") %>% pull(v_identity) -> x
v_identity %>% subset(type == "VH") %>% pull(v_identity) -> y
t_test <- t.test(x, y)
p_identity <- ggplot(v_identity, aes(x =  type, y = v_identity)) +
  geom_violin(aes(fill = type)) +
  geom_boxplot(width=0.1) + 
  # geom_jitter(width = 0.2) +
  theme_classic() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 0, colour = "black", hjust=0.95,vjust=0.3),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "", y = "Identity", color = "Type",  caption = paste0("p_value = ", t_test$p.value))
p_identity 


pdf("results/v_cdna/v_alignment.pdf")
p_v_subfamily
p_identity 
dev.off()
