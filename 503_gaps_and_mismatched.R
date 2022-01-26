

load("data/align.Rda")
load("data/v_cdna/v_cdna.Rda")
load("data/v_gdna/v_gdna.Rda")
source("R/find_hotspots.R")


v_germline <- v_gdna %>% select(id, imgt_aa, imgt_nt, type) %>% 
  rename(v_call = id, imgt_aa_germline = imgt_aa, imgt_nt_germline = imgt_nt, type_germline = type)


v_pairs <- v_cdna %>% 
  separate(v_call, "v_call", sep = ",", extra = "drop") %>%
  select(sequence_id, v_call, imgt_aa, imgt_nt, type) %>%
  rename(imgt_aa_rearranged = imgt_aa, imgt_nt_rearranged = imgt_nt, type_rearranged = type) %>%
  left_join(v_germline) %>%
  subset(!is.na(imgt_aa_germline) & !is.na(imgt_aa_rearranged)) %>%
  subset(type_rearranged == type_germline)

total_count <- nrow(v_pairs)

vs=v_pairs
vhhs <- vs %>% subset(type_germline == "VHH") 
vhs <- vs %>% subset(type_germline == "VH") 




# Gaps and inserts
indels_vhh <- map2(vhhs$imgt_aa_germline, vhhs$imgt_aa_rearranged, detect_gaps) %>%
  bind_rows(.id = "no") %>% mutate(type = "VHH")
indels_vh <- map2(vhs$imgt_aa_germline, vhs$imgt_aa_rearranged, detect_gaps) %>%
  bind_rows(.id = "no") %>% mutate(type = "VH")

# Mismatches
mismatches_vhh <- map2(vhhs$imgt_aa_germline, vhhs$imgt_aa_rearranged, detect_mismatches) %>%
  bind_rows(.id = "no") %>% mutate(variant = "mismatch") %>% mutate(type = "VHH")
mismatches_vh <- map2(vhs$imgt_aa_germline, vhs$imgt_aa_rearranged, detect_mismatches) %>%
  bind_rows(.id = "no") %>% mutate(variant = "mismatch") %>% mutate(type = "VH")



variants <- bind_rows(indels_vhh, indels_vh, mismatches_vhh, mismatches_vh) %>%
  group_by(type, variant, site) %>%
  summarise(count = n()) %>%
  mutate(region = case_when(site < 27 ~ "FR1",
                            site < 39 ~ "CDR1",
                            site < 56 ~ "FR2",
                            site < 66 ~ "CDR2",
                            site < 105 ~ "FR3")) %>%
  mutate(count = ifelse(type == "VH", -count, count)) %>%
  mutate(prop = count/total_count)


p_mismatches <- variants %>% subset(variant == "mismatch") %>%
  # ggplot(aes(x = site, y = prop, fill = region)) +
  ggplot(aes(x = site, y = prop, fill = as.factor(site))) +
  geom_hline(yintercept = 0) +
  geom_bar(stat = "identity", width = 0.6) +
  facet_grid(variant ~ ., scales = "free", space = "free", shrink = FALSE) +
  scale_x_continuous(breaks = seq(0, 105, 1)) +
  scale_y_continuous(limits = c(-0.4, 0.4)) +
  theme_classic() +
  # theme(legend.position = "top", strip.background = element_blank()) +
  # labs(x = "", y = "", fill = "region", caption =  "VHH/VH") +
  theme(legend.position = "none", strip.background = element_blank()) +
  labs(x = "", y = "", caption =  "VHH/VH") 
p_mismatches

p_indels <- variants %>% subset(variant != "mismatch") %>%
  # ggplot(aes(x = site, y = prop, fill = region), width = 0.6) +
  ggplot(aes(x = site, y = prop, fill = as.factor(site)), width = 0.6) +
  geom_hline(yintercept = 0) +
  geom_bar(stat = "identity",width = 0.6) +
  facet_grid(variant ~ ., scales = "free", space = "free", shrink = FALSE) +
  scale_y_continuous(limits = c(-0.06, 0.06)) +
  scale_x_continuous(limits = c(0, 105), breaks = seq(0, 105, 1)) +
  theme_classic() +
  # theme(legend.position = "top", strip.background = element_blank()) +
  # labs(x = "", y = "", fill = "region", caption =  "VHH/VH") +
  theme(legend.position = "nonep", strip.background = element_blank()) +
  labs(x = "", y = "", caption =  "VHH/VH")
p_indels 


pdf("results/v_cdna/mutation_hotspots.pdf", width = 6, height = 4)
p_mismatches
p_indels 
dev.off()



