#!/usr/bin/env Rscript

# Load DJC sequences



load("data/djc_cdna/vdjc.Rda")



# ## CDR3 AA composite
# max_length <- max(nchar(vdj$cdr3_aa), na.rm = TRUE) + 1
# aa_composite <- vdjc_d %>% select(species, sample, cdr3_aa, count) %>%
#   subset(!is.na(cdr3_aa)) %>%
#   separate(cdr3_aa, into = paste0("a", c(1:max_length)), sep = "") %>%
#   gather(key = "aa", value = "composite", -c(species, sample, count)) 
# 
# d <- aa_composite %>% subset(composite != "") %>%
#   group_by(species, sample, composite) %>% summarise(count = sum(count)) %>%
#   group_by(species, sample) %>% add_tally(count, name = "total_count") %>%
#   mutate(proportion = count / total_count)
# 
# p_aa_composite <- ggplot(d, aes(x = sample, y = proportion, fill = composite)) +
#   geom_bar(stat = "identity") +
#   theme_classic() +
#   labs(x = "", y = "Proportion", fill = "Composite")
# p_aa_composite
# 
# pdf("results/figure_4_cdr3_aa_composite.pdf",  width = 5, height = 4)
# p_aa_composite
# dev.off()




##
vdjc %>%
  subset(!is.na(fwr4_aa)) %>% subset(!is.na(fwr3_aa)) %>%
  select(sample, species, cdr3) %>% 
  mutate(cystein = str_count(cdr3, "C")) %>%
  group_by(sample, species, cystein) %>%
  summarise(count = n()) %>% 
  drop_na() %>%
  group_by(sample) %>%
  add_tally(count, name = "total_count") %>%
  mutate(prop = count / total_count)


n_cys <- vdjc %>%
  subset(!is.na(fwr4_aa)) %>% subset(!is.na(fwr3_aa)) %>%
  select(cdr3, fwr4_aa) %>%
  rownames_to_column("id")%>%
  gather("region", "seq", -id) %>%
  mutate(length = nchar(seq)) %>%
  mutate(n_cys = str_count(seq, "C")) %>%
  group_by(region, n_cys) %>% summarise(count = n()) %>%
  group_by(region) %>% mutate(prop = count / sum(count)) 

p_cdr3_cys <- n_cys %>%
  mutate(label = paste0(n_cys, " (", round(prop *100, 2), "%)")) %>%
  mutate(n_cys = as.factor(n_cys)) %>%
  subset(region == "cdr3") %>%
  ggplot(aes(x = 1, y = prop, fill = reorder(label, -prop))) +
  geom_bar(stat = "identity") +
  # xlim(0, 2) +
  coord_polar(theta = "y") +
  theme_void()


p_fwr4_cys <- n_cys %>%
  mutate(label = paste0(n_cys, " (", round(prop *100, 2), "%)")) %>%
  mutate(n_cys = as.factor(n_cys)) %>%
  subset(region == "fwr4_aa") %>%
  ggplot(aes(x = 1, y = prop, fill = reorder(label, -prop))) +
  geom_bar(stat = "identity") +
  # xlim(0, 2) +
  coord_polar(theta = "y") +
  theme_void()



pdf("results/djc_cdna/cysteins_2.pdf")
p_cdr3_cys 
p_fwr4_cys 
dev.off()

# FR4
fr4 <- vdjc %>% select(sample, species, fr4_seq) %>% 
  mutate(cystein = str_count(fr4_seq, "C")) %>%
  group_by(sample, species, cystein) %>%
  summarise(count = n()) %>% 
  drop_na() %>%
  group_by(sample) %>%
  add_tally(count, name = "total_count") %>%
  mutate(prop = count / total_count)

p_cystein_fr4 <- ggplot(fr4, aes(x = sample, y = prop, fill = as.factor(cystein))) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, colour = "black", hjust = 0.95, vjust = 0.3),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "", y = "Length", fill = "Species")

p_cystein_fr4



pdf("results/djc_cdna/cysteins.pdf")
p_cystein_fr4
p_cystein_cdr3
dev.off()


## Chi-square-test

