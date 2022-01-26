#!/usr/bin/env Rscript
#


load("data/v_cdna/v_cdna.Rda")



full_length <- v_cdna %>%
  mutate(full = nchar(str_remove_all(imgt_aa, "\\."))) %>%
  group_by(species, sample, type, full) %>% summarise(count = n()) %>%
  group_by(species, sample, type) %>% mutate(prop = count / sum(count))

p_full_length <- full_length %>%
  ggplot(aes(x = full , fill = sample)) +
  geom_density(alpha = 0.3) +
  facet_grid(type ~ .) +
  theme_classic()

pdf("results/v_cdna/full_length.pdf", width = 8, height = 5)
p_full_length
dev.off()

#### Region length
region_length <- v_cdna %>%
  mutate(
    fr1 = nchar(fwr1_aa),
    fr2 = nchar(fwr2_aa),
    fr3 = nchar(fwr3_aa),
    cdr1 = nchar(cdr1_aa),
    cdr2 = nchar(cdr2_aa)
  ) %>%
  select(species, sample, type, fr1, fr2, fr3, cdr1, cdr2) %>%
  gather(key = "region", value = "length", -c(species, sample, type)) 

p_region_region <- region_length %>%
  group_by(sample, type, region, length) %>% summarise(count = n()) %>%
  group_by(sample, type, region) %>% mutate(prop = count / sum(count)) %>%
  mutate(length = as.factor(length)) %>%
  mutate(region = factor(region, levels = c("fr1", "cdr1", "fr2", "cdr2", "fr3"))) %>%
  ggplot(aes(x = sample, y = prop, group = prop, fill = fct_reorder(length, prop, .fun = sum))) +
  geom_bar(stat = "identity") +
  facet_grid(type ~ region) +
  theme_classic() +
  theme(strip.background = element_blank(),
        # strip.text.x = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 0, hjust = 0.95, vjust = 0.2),
        legend.position = "right") +
  labs(x = "", y = "Proportion", fill = "Length") +
  guides(fill=guide_legend(ncol = 2, byrow = TRUE))
p_region_region

pdf("results/v_cdna/region_length_by_sample_and_type.pdf", width = 8, height = 5)
p_region_region
dev.off()



region_length %>%
  group_by(sample, type, region, length) %>% summarise(count = n()) %>%
  group_by(sample, type, region) %>% mutate(prop = count / sum(count)) %>%
  mutate(length = as.factor(length)) %>%
  subset(region %in% c("cdr1", "cdr2")) %>%
  # mutate(region = factor(region, levels = c("fr1", "cdr1", "fr2", "cdr2", "fr3"))) %>%
  ggplot(aes(x = sample, y = prop, group = prop, fill = fct_reorder(length, prop, .fun = sum))) +
  geom_bar(stat = "identity") +
  facet_grid(type ~ region) +
  theme_classic() +
  theme(strip.background = element_blank(),
        # strip.text.x = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 0, hjust = 0.95, vjust = 0.2),
        legend.position = "right") +
  labs(x = "", y = "Proportion", fill = "Length") +
  guides(fill=guide_legend(ncol = 1, byrow = TRUE))



region_length %>%
  group_by(type, region, length) %>% summarise(count = n()) %>%
  group_by(type, region) %>% mutate(prop = count / sum(count)) %>%
  mutate(length = as.factor(length)) %>%
  subset(region %in% c("cdr1", "cdr2")) %>%
  # mutate(region = factor(region, levels = c("fr1", "cdr1", "fr2", "cdr2", "fr3"))) %>%
  ggplot(aes(x = type, y = prop, group = prop, fill = fct_reorder(length, prop, .fun = sum))) +
  geom_bar(stat = "identity") +
  facet_grid(type ~ region) +
  theme_classic() +
  theme(strip.background = element_blank(),
        # strip.text.x = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 0, hjust = 0.95, vjust = 0.2),
        legend.position = "right") +
  labs(x = "", y = "Proportion", fill = "Length") +
  guides(fill=guide_legend(ncol = 1, byrow = TRUE))


region_length %>%
  subset(region %in% c("cdr1", "cdr2")) %>%
  nest_by(type, region) %>%
  mutate(model = list(chisq.test(data$sample, data$length))) %>%
  summarise(broom::tidy(model))




region_length %>%
  subset(region %in% c("cdr1", "cdr2")) %>%
  subset(type == "VH") %>% subset(region  == "cdr1") %>%
  select(-type, -region, -species) %>% 
  table() %>% chisq.test()

region_length %>%
  subset(region %in% c("cdr1", "cdr2")) %>%
  group_by(type, region) %>% summarise(avg = mean(length), sd = sd(length)) %>%
  ggplot(aes(x = type, y = avg, fill = type)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = avg-sd, ymax = avg+sd)) +
  facet_grid(.~region) +
  theme_classic()




# plot_region_length <- function(subregion) {
#   region_length %>% 
#     group_by(sample, type, region, length) %>% summarise(count = n()) %>%
#     group_by(sample, type, region) %>% mutate(prop = count / sum(count)) %>%
#     subset(region == subregion) %>%
#     mutate(length = as.factor(length)) %>%
#     ggplot(aes(x = type, y = prop, fill = length)) +
#     geom_bar(stat = "identity") +
#     facet_grid(. ~ sample, scales = "free", space = "free") + 
#     theme_classic() +
#     theme(strip.background = element_blank(),
#           # strip.text.x = element_text(face = "bold", size = 12),
#           axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2),
#           legend.position = "right") +
#     labs(x = "", y = "Proportion") +
#     guides(fill=guide_legend(ncol = 2, byrow = TRUE))
# }
# 
# p_cdr1 <- plot_region_length("cdr1")
# p_cdr2 <- plot_region_length("cdr2")
# p_fr1 <- plot_region_length("fr1")
# p_fr2 <- plot_region_length("fr2")
# p_fr3 <- plot_region_length("fr3")
# 
# gridExtra::grid.arrange(p_cdr1, p_cdr2, p_fr1, p_fr2, p_fr3, ncol = 1)
# 
# pdf("results/v_cdna/region_length_by_sample_and_type_2.pdf", width = 6, height = 12)
# gridExtra::grid.arrange(p_cdr1, p_cdr2, p_fr1, p_fr2, p_fr3, ncol = 1)
# dev.off()





#### Compare region length by type
p_region_lenth_by_type <- region_length %>% subset(region %in% c("cdr1", "cdr2")) %>%
  group_by(type, region, length) %>% summarise(count = n()) %>%
  group_by(type, region) %>% mutate(prop = count / sum(count)) %>%
  mutate(length = as.factor(length)) %>%
  mutate(region = factor(region, levels = c("fr1", "cdr1", "fr2", "cdr2", "fr3"))) %>%
  ggplot(aes(x = type, y = prop, group = prop, fill = fct_reorder(length, prop, .fun = sum))) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ region) +
  theme_classic() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.95, vjust = 0.2),
        legend.position = "right") +
  labs(x = "", y = "Proportion") +
  guides(fill=guide_legend(ncol = 2, byrow = TRUE,title = "Length"))
p_region_lenth_by_type




region_length %>% subset(region %in% c("cdr1", "cdr2")) %>%
  nest_by(region) %>%
  mutate(model = list(chisq.test(data$type, data$length))) %>%
  summarise(broom::tidy(model))

chisq_test <- region_length %>% subset(region %in% c("cdr1", "cdr2")) %>%
  group_by(type, region, length) %>% summarise(count = n()) %>%
  subset(region == "cdr1") %>% ungroup() %>% select(-region) %>%
  spread(length, count, fill = 0) %>%
  column_to_rownames("type") %>% chisq.test()

# p-value < 2.2e-16

#### Compare region length by species and type
## Plot
p_region_lenth_by_species_and_type <- region_length %>% subset(region %in% c("cdr1", "cdr2")) %>%
  group_by(species, type, region, length) %>% summarise(count = n()) %>%
  group_by(species, type, region) %>% mutate(prop = count / sum(count)) %>%
  mutate(length = as.factor(length)) %>%
  mutate(region = factor(region, levels = c("fr1", "cdr1", "fr2", "cdr2", "fr3"))) %>%
  ggplot(aes(x = species, y = prop, group = prop, fill = fct_reorder(length, prop, .fun = sum))) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ type + region) +
  theme_classic() +
  theme(strip.background = element_blank(),
        # strip.text.x = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2),
        legend.position = "right") +
  labs(x = "", y = "Proportion") +
  guides(fill=guide_legend(ncol = 2, byrow = TRUE,title = "Length"))
p_region_lenth_by_species_and_type

## Compare by species
region_by_species_and_type


v_cdna %>%
  mutate(
    fr1 = nchar(fwr1_aa),
    fr2 = nchar(fwr2_aa),
    fr3 = nchar(fwr3_aa),
    cdr1 = nchar(cdr1_aa),
    cdr2 = nchar(cdr2_aa)
  ) %>%
  select(species, type, cdr1, cdr2) %>%
  gather(key = "region", value = "length", -c(species, type)) %>%
  nest_by(region, type) %>%
  mutate(model = list(chisq.test(data$species, data$length))) %>%
  summarise(broom::tidy(model)) %>%
  set_names(str_replace_all(names(.), "\\.", "_")) %>%
  rename(group = type) %>%
  mutate(by = "dome_vs_wild")


v_cdna %>%
  mutate(
    fr1 = nchar(fwr1_aa),
    fr2 = nchar(fwr2_aa),
    fr3 = nchar(fwr3_aa),
    cdr1 = nchar(cdr1_aa),
    cdr2 = nchar(cdr2_aa)
  ) %>% select(species, type, cdr1, cdr2) %>%
  gather(key = "region", value = "length", -c(species, type)) %>%
  subset(type == "VHH") %>% subset(region == "cdr1") %>%
  select(species, length) %>% 
  table() %>% chisq.test()




# chisq_test <- bind_rows(chisq_test_by_species, chisq_test_by_type)

# Save
pdf("results/v_cdna/region_length_by_species_and_type.pdf")
p_region_lenth_by_type
p_region_lenth_by_species_and_type
dev.off()
write.csv(chisq_test, "results/v_cdna/region_length_chisq_test.csv")

