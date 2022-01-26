
load("data/djc_cdna/vdjc.Rda")

v_reference_extend <- read_csv("imported/refs/v_reference_extend.csv")

cdr3_camelids <-  v_reference_extend %>% rename_all(tolower) %>%
  mutate(type = type_adjusted) %>%
  select(id, from, species, molecule,type, cdr3) %>%
  mutate(subgene = "IGHV3", family = "camelids") %>% 
  subset(molecule=="cDNA") %>%
  select(species, type, cdr3) %>%
  subset(species == "Lama glama") %>%
  subset(!is.na(cdr3)) %>%
  mutate(seq = str_remove_all(cdr3, "\\.")) %>%
  mutate(cdr3_length = nchar(seq)) %>%
  select(species, type, cdr3_length)


cdr3_bactrian <- vdjc %>% select(species,  cdr3_length) %>%
  mutate(type = "VHH", species = "Camelus bactrianus") 

  
p_length_density <- bind_rows(cdr3_bactrian, cdr3_camelids) %>%
  ggplot(aes(x = cdr3_length, group = species)) +
  # geom_bar(stat = "identity", aes(fill = sample)) +
  geom_density(aes(y = ..density.., fill = species), lwd = 1, alpha = 0.5, adjust = 2) +
  facet_grid(type ~ .) +
  theme_classic() +
  labs(fill = "Sample", color = "Sample", x = "Length", y = "Porportion")
p_length_density 

# d <- vdjc %>% 
#   group_by(sample, species, cdr3_length) %>% summarise(count = n()) %>%
#   group_by(sample) %>% mutate(prop = count / sum(count)) %>%
#   arrange(sample, -prop)


# # Density
# p_cdr3_length_density <- ggplot(d, aes(x = cdr3_length, y = prop, group = sample)) +
#   # geom_bar(stat = "identity", aes(fill = sample)) +
#   geom_density(aes(y = ..density.., color = sample), lwd = 1) +
#   # facet_grid(sample ~ .) +
#   theme_classic() +
#   labs(fill = "Sample", color = "Sample", x = "Length", y = "Porportion")
# p_cdr3_length_density 
# 




p_bar <- bind_rows(cdr3_bactrian, cdr3_camelids) %>%
  group_by(species, type) %>% 
  summarise(mean = mean(cdr3_length), sd = sd(cdr3_length)) %>%
  ggplot(aes(x = type, y = mean, fill = type)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(x=type, ymin=mean-sd, ymax=mean+sd), width = .5) +
  facet_grid(.~species)+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.95, vjust = 0.2))
p_bar

pdf("results/djc_cdna/cdr3_length_2.pdf")
p_length_density 
p_bar
dev.off()



cdr3_camelids <-  v_reference_extend %>% rename_all(tolower) %>%
  mutate(type = type_adjusted) %>%
  select(id, from, species, molecule,type, cdr3) %>%
  mutate(subgene = "IGHV3", family = "camelids") %>% 
  subset(molecule=="cDNA") %>%
  select(species, type, cdr3) %>%
  subset(species != "Camelus dromedarius") %>%
  subset(!is.na(cdr3))

cdr3_bactrian <- vdjc %>% select(species,  cdr3) %>%
  mutate(type = "VHH", species = "Camelus bactrianus") 



n_cys <- bind_rows(cdr3_bactrian, cdr3_camelids) %>%
  subset(nchar(cdr3)>0) %>%
  mutate(n_cys = str_count(cdr3, "C")) %>%
  mutate(group = ifelse(n_cys > 0, ">=1", "0")) %>%
  group_by(species, type, group) %>% summarise(count = n()) %>%
  group_by(species, type) %>% mutate(prop = count/sum(count)) 


p_n_cys <- n_cys %>%
  mutate(type = factor(type, levels = c("VHH", "VH"))) %>%
  # subset(subgene == "IGHV3") %>%
  # ggplot(aes(x = species, y = prop, fill = species)) +
  ggplot(aes(x = type, y = prop)) +
  geom_bar(stat = "identity", aes(fill = group)) +
  facet_grid(. ~ species, scales = "free",space = "free") +
  theme_classic() 
  theme(axis.text.x = element_text(angle = 90, 
                                   colour = "black", 
                                   hjust = 0.95, 
                                   vjust = 0.2)) +
  labs(x = "", y = "Proportion", fill = "Species") 

p_n_cys


pdf("results/djc_cdna/cdr3_cystein.pdf")
p_n_cys
dev.off()



# Boxplot
p_cdr3_length_boxplot <- ggplot(d, aes(x = sample, y = cdr3_length, fill = species)) +
  geom_boxplot() +
  theme_classic() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, colour = "black", hjust = 0.95, vjust = 0.3),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "", y = "Length", fill = "Species")
 
p_cdr3_length_boxplot

pdf("results/djc_cdna/cdr3_length.pdf")
p_cdr3_length_density 
p_cdr3_length_boxplot
dev.off()








v_reference_extend <- read_csv("imported/refs/v_reference_extend.csv")

v_camelids <-  v_reference_extend %>% rename_all(tolower) %>%
  mutate(type = type_adjusted) %>%
  select(id, from, species, molecule,type, cdr3) %>%
  mutate(subgene = "IGHV3", family = "camelids") 


d <- v_camelids %>% 
  subset(molecule=="cDNA") %>%
  select(species, type, cdr3) %>%
  subset(species != "Camelus dromedarius") %>%
  mutate(seq = str_remove_all(cdr3, "\\.")) %>%
  mutate(length = nchar(seq)) 

ggplot(d, aes(x = length, group = species)) +
  # geom_bar(stat = "identity", aes(fill = sample)) +
  geom_density(aes(y = ..density.., fill = species), lwd = 1, alpha = 0.5) +
  facet_grid(type ~ .) +
  theme_classic() +
  labs(fill = "Sample", color = "Sample", x = "Length", y = "Porportion")

