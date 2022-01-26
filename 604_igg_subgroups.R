# IgG subgroups


source("R/analyse_subgroups.R")
load("data/djc_cdna/vdjc.Rda")

# library(ape)
# library(ggmsa)
# library(Biostrings)




# Check hinge sequence and determine subgroups
hinge <- vdjc %>% 
  subset(stop_codon == FALSE) %>%
  subset(!is.na(h_seq)) %>%
  group_by(h_seq) %>% summarise(count = n()) %>%
  mutate(h_seq = str_remove_all(h_seq, "-")) %>%
  arrange(-count) %>%
  rownames_to_column("seq_id") %>%
  mutate(seq_id = paste0("h_", seq_id)) 

subgroups <- hinge %>% head(5) %>% select(h_seq) %>%
  mutate(subgroup = c("IgG2a", "IgG3-gg", "IgG3-ev", "IgG2c", "IgG1a"))
save(subgroups, file = "data/djc_cdna/subgroups.Rda")



# Summary subgroups
subgroups_individuals <- vdjc %>% 
  left_join(subgroups, by = "h_seq") %>% subset(!is.na(subgroup)) %>%
  subset(subgroup != "IgG1a") %>%
  # separate(subgroup, into=c("subgroup", NA), remove = F, sep="-") %>%
  group_by(sample, subgroup) %>% summarise(count = n()) %>%
  group_by(sample) %>% mutate(prop = count / sum(count)) %>%
  arrange(sample, -prop) %>%
  mutate(label = paste0(round(prop *100 , 2), "%")) 
  
# Plot heatmap
subgroup_heatmap(subgroups_individuals, "results/djc_cdna/heatmap_hinge_subgroups.pdf")






subgroups_individuals %>%
  group_by(subgroup) %>% summarise(prop = mean(prop)) %>%
  mutate(prop / sum(prop))


subgroups_individuals %>%
  group_by(subgroup) %>% summarise(prop = mean(prop))


# Overall subgroup
# Pie plot
subgroups_overall <- subgroups_individuals %>%
  group_by(subgroup) %>% summarise(count = sum(count)) %>%
  mutate(prop = round(count / sum(count) * 100, 4))

pdf("results/djc_cdna/pie_subgroup.pdf")
subgroup_pie(subgroups_overall)
dev.off()




subgroups_individuals %>%
  ggplot(aes(x=sample, y = prop, fill = fct_reorder(subgroup, prop, .fun = sum))) +
  geom_bar(stat = "identity") +
  theme_classic() 
  facet_grid(.~sample, scales = "free")


# Compare by wild and domestic
# Chi square test 
chisq_test <- subgroups_individuals %>% 
  mutate(species = ifelse(sample %in% c("cDNA-wild-1", "cDNA-wild3"), "wild", "domestic")) %>%
  group_by(species, subgroup) %>% summarise(count = sum(count)) %>%
  spread(species, count) %>%
  column_to_rownames("subgroup") %>% chisq.test()

d <- subgroups_individuals %>% 
  mutate(species = ifelse(sample %in% c("cDNA-wild-1", "cDNA-wild3"), "wild", "domestic")) %>%
  group_by(species, subgroup) %>% summarise(count = sum(count)) %>%
  group_by(species) %>% mutate(prop = count / sum(count)) %>%
  mutate(subgroup = factor(subgroup)) %>%
  mutate(subgroup = fct_reorder(subgroup, -prop, .fun=sum)) %>%
  mutate(lab = paste0(round(prop, 4) *100, "%"))

subgroup_by_specis <- ggplot(d, aes(x = species, y = prop, fill = subgroup)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = lab), position = position_fill(vjust = 0.5, reverse = F)) +
  labs(caption = as.character(chisq_test$p.value)) +
  theme_classic()
subgroup_by_specis
pdf("results/djc_cdna/bar_subgroup_by_species.pdf")
subgroup_by_specis
dev.off()

