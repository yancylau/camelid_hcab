library(ggalluvial)




source("R/colors.R")
load("data/djc_cdna/vdjc.Rda")
load("data/djc_cdna/subgroups.Rda")




d <- vdjc %>% left_join(subgroups, by = "h_seq") %>%
  # left_join(ch2_subgroups, by = "ch2_seq") %>%
  subset(!is.na(subgroup) & !is.na(subgroup)) %>%
  select(sample, d_call, j_call, subgroup) %>%
  group_by_all() %>% summarise(count = n()) %>%
  subset(subgroup != "IgG1a") %>%
  subset(!is.na(d_call))


p_interaction <- d %>%
  mutate(d_call = recode(d_call, 
                         IGHD2 = "d1", 
                         IGHD1 = "d2", 
                         IGHD5 = "d3", 
                         IGHD3 = "d4", 
                         IGHD4 = "d5", 
                         IGHD6 = "d6", 
                         IGHD7 = "d7")) %>%
  mutate(j_call = recode(j_call, 
                         IGHJ4 = "j1", 
                         IGHJ6= "j2", 
                         IGHJ7 = "j3", 
                         IGHJ3 = "j4", 
                         IGHJ2 = "j5", 
                         IGHJ5 = "j6", 
                         `IGHJ1/ORF` = "j7")) %>%
  mutate(subgroup = recode(subgroup,
                           IgG2a = "s1",
                           `IgG3-gg` = "s2",
                           `IgG3-ev` = "s3",
                           IgG2c = "s4")) %>%
  ggplot(aes(axis1= d_call, axis2 = j_call, axis3 = subgroup, y = count)) +
  scale_x_discrete(limits = c("D", "J,", "Subgroup"), expand = c(.2, .05))+
  geom_alluvium(aes(fill = sample), alpha = 0.5) +
  geom_stratum()  +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  # scale_fill_manual(values = colors) +
  theme_void() +
  labs(x = "", y = "Count", fill = "Sample")
  
p_interaction


pdf("results/djc_cdna/interaction.pdf")
p_interaction
dev.off()


