#!/usr/bin/env Rscript
#
#
# Determine sequence type based on hallmarks
# Site      VH      VHH
#  42     V,I      F,Y
#  49       G      E,Q
#  50       L      C,R
#  52     F,W      G,L
#
# Reference:
# Muyldermans, S. (2013).  Annu. Rev. Biochem. 82, 775–797.


# iglast: 56,011
# Determined by hallmarks: 44542 (79.52366%), undetermined: 11469 (20.47634%)
#
# VH: 30417 (68.28836%) + VHH: 14125 (31.71164%) 
#
#            mix, novel traditional unknown
# unknown     0,     0           0     186
# VH        914,  3527       30417       0
# VHH      2148,  4694       14125       0

load("data/v_cdna/igblast.Rda")


## Determine sequences type based on hallmarks
v_cdna_hallmarks <- igblast %>%
  mutate(a42 = substr(imgt_aa, 42, 42), 
         a49 = substr(imgt_aa, 49, 49), 
         a50 = substr(imgt_aa, 50, 50), 
         a52 = substr(imgt_aa, 52, 52)) %>% 
  mutate(n42 = case_when(a42 %in% c("V", "I") ~ 1,
                         a42 %in% c("F", "Y") ~ -1,
                         TRUE ~ 0),
         n49 = case_when(a49 %in% c("G") ~ 1,
                         a49 %in% c("E", "Q") ~ -1,
                         TRUE ~ 0),
         n50 = case_when(a50 %in% c("L") ~ 1,
                         a50 %in% c("C", "R") ~ -1,
                         TRUE ~ 0),
         n52 = case_when(a52 %in% c("F", "W") ~ 1,
                         a52 %in% c("G", "L") ~ -1,
                         TRUE ~ 0)) %>%
  mutate(sum = n42 + n49 + n50 + n52, abs_sum = abs(n42) + abs(n49) + abs(n50) + abs(n52)) %>%
  mutate(type1 = ifelse(sum > 0, "VH", ifelse(sum < 0, "VHH", "unknown"))) %>%
  mutate(type2 = case_when(sum == 4 ~ "traditional",                          # VH-traditional
                           sum == 3 ~ "novel",                                # VH-novel
                           sum == 2 & abs_sum == 2 ~ "novel",                 # VH-novel
                           sum == 2 & abs_sum == 4 ~ "mix",                   # VH-mix
                           sum == 1 & abs_sum == 1 ~ "novel",                 # VH-novel
                           sum == 1 & abs_sum == 3 ~ "mix",                   # VH-mix
                           sum == -4 ~ "traditional",                         # VHH-traditional
                           sum == -3 ~ "novel",                               # VHH-novel
                           sum == -2 & abs_sum == 2 ~ "novel",                # VHH-novel
                           sum == -2 & abs_sum == 4 ~ "mix",                  # VHH-mix
                           sum == -1 & abs_sum == 1 ~ "novel",                # VHH-novel
                           sum == -1 & abs_sum == 3 ~ "mix",                  # VHH-mix
                           TRUE ~ "unknown")) %>%
  #mutate(type = factor(type1, levels = c("VH", "VHH", "unknown"))) %>%              # Set types as factor and its levels
  mutate(subtype = factor(type2, levels = c("traditional", "novel", "mix", "unknown"))) %>%
  mutate(type = ifelse(type2 == "traditional", type1, "unknown"))


save(v_cdna_hallmarks, file = "data/v_cdna/v_cdna_hallmarks.Rda")


## Summary sequences types
# Canonical: VH/VH
d1 <- v_cdna_hallmarks %>% 
  group_by(species, sample, type) %>% summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  mutate(label = paste0(round(proportion, 4) * 100, "%")) %>%
  mutate(level = "type") %>%
  rename(group = type)

# non-canonical: traditional, mix, noval
d2 <- v_cdna_hallmarks %>% 
  group_by(species, sample, type, subtype) %>% summarise(count = n()) %>%
  group_by(sample) %>% mutate(proportion = count / sum(count)) %>%
  mutate(label = paste0(round(proportion, 4) * 100, "%")) %>%
  mutate(group = paste(type, subtype)) %>%
  mutate(level = "subtype") %>%
  select(species, sample, count, proportion, label, group, level) %>%
  mutate(group = factor(group, levels = c("VH novel", "VH mix", "VH traditional", 
                                          "VHH traditional", "VHH mix", "VHH novel", 
                                          "unknown unknown")))

p_types_sample <- ggplot() + 
  geom_col(data = d1, aes(x = 3, y = proportion, fill = group), color = "white", size = 0) +
  geom_text(data = d1, aes(x = 3, y = proportion, fill = group, label = label), position = position_stack(vjust = 0.5), size = 3) +
  labs(fill = "Type") +
  ggnewscale::new_scale_fill() +
  geom_col(data = d2, aes(x = 2, y = proportion, fill = group), color = "white", size = 0) +
  geom_text(data = d2, aes(x = 2, y = proportion, fill = group, label = label), position = position_stack(vjust = 0.5), size = 2) +
  xlim(0, 3.5) + labs(x = NULL, y = NULL) +
  coord_polar(theta = "y") +
  theme_void() +
  facet_wrap(. ~ sample, nrow = 2) +
  labs(fill = "Subtype")
p_types_sample


d1 <- v_cdna_hallmarks %>% 
  group_by(type) %>% summarise(count = n()) %>%
  ungroup() %>% mutate(proportion = count / sum(count)) %>%
  mutate(label = paste0(type, " (", round(proportion, 4) * 100, "%)")) %>%
  mutate(level = "type") %>%
  rename(group = type) 

# non-canonical: traditional, mix, novel
d2 <- v_cdna_hallmarks %>% 
  group_by(type, subtype) %>% summarise(count = n()) %>%
  ungroup() %>% mutate(proportion = count / sum(count)) %>%
  mutate(group = paste(type, subtype)) %>%
  mutate(label = paste0(group, " (", round(proportion, 4) * 100, "%)")) %>%
  mutate(level = "subtype") %>%
  select(count, proportion, label, group, level) %>%
  mutate(group = factor(group, levels = c("VH novel", "VH mix", "VH traditional", 
                                          "VHH traditional", "VHH mix", "VHH novel", 
                                          "unknown unknown")))

p_types <- ggplot() + 
  geom_col(data = d1, aes(x = 3, y = proportion, fill = label), color = "white", size = 0) +
  #geom_text(data = d1, aes(x = 3, y = proportion, fill = label, label = label), position = position_stack(vjust = 0.5), size = 3) +
  labs(fill = "Type") +
  ggnewscale::new_scale_fill() +
  geom_col(data = d2, aes(x = 2, y = proportion, fill = label), color = "white", size = 0) +
  #geom_text(data = d2, aes(x = 2, y = proportion, fill = label, label = label), position = position_stack(vjust = 0.5), size = 2) +
  xlim(0, 3.5) + labs(x = NULL, y = NULL) +
  coord_polar(theta = "y") +
  theme_void() 

p_types







