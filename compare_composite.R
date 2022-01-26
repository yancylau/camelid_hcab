

# Calculate aa composite on FR sequences



## Calculate aa frequency  -----------------------------------------------------
calculate_allele_frequency <- function(pos) {
  frs %>%
    mutate(composite = substr(frs, pos, pos)) %>%
    select(type, composite) %>%
    subset(composite != "X") %>%                     # Remove trancated sequnces
    # mutate(composite = ifelse(composite %in% c(".", "-", "X"), NA, composite)) %>%
    group_by(type, composite) %>% summarise(count = n()) %>%
    group_by(type) %>% add_tally(count, name = "total_count") %>%
    mutate(prop = count / total_count) %>%
    mutate(pos = pos) %>%
    select(pos, composite, type, count, prop)
}



site_dict <- tibble(pos = c(1:82), site = c(1:26, 39:55, 66:104))
allele_frequency <- map(1:82, calculate_allele_frequency) %>% bind_rows() %>%
  left_join(site_dict, by = "pos") %>%
  mutate(region = case_when(site < 27 ~ "FR1", 
                            site < 56 ~ "FR2", 
                            site < 105 ~ "FR3")) %>%
  identity()





## Plot site -------------------------------------------------------------------
p_allele_frequency <- allele_frequency %>%
  mutate(site = factor(site)) %>%
  mutate(prop = ifelse(type == "VHH", prop, -prop)) %>%
  mutate(type = factor(type, levels = c("VHH", "VH"))) %>%
  # ggplot(aes(x = site, y = prop, fill = composite)) +
  ggplot(aes(x = site, y = prop, group = count, fill = fct_reorder(composite, prop, .fun = sum))) +
  # geom_segment(aes(xend = composite, yend = 0, color = composite), size = 1) +
  # geom_point(aes(color = composite ), size = 4) +
  geom_hline(yintercept = 0, color = "darkgrey") +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  theme_classic() +
  theme(strip.background = element_blank(),
        # strip.text.x = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2),
        legend.position = "bottom") +
  facet_grid(. ~ region, scales = "free_x", space = "free") +
  labs(x = "", y = "Proportion") +
  guides(fill=guide_legend(nrow = 2, byrow = TRUE))
p_allele_frequency


