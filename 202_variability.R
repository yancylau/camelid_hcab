#!/usr/bin/env Rscript

# Calculate variability and compare all sites


dir.create("results/camelids", showWarnings = F)


load("data/v_gdna/v_gdna.Rda")
load("data/v_gdna/otus.Rda")
load("data/refs/v_camelids.Rda")



## Function: Calculate variability index for one site
calculate_variability <- function(seq) {
  # Remove gaps or in this site
  seq <- seq[! seq %in% c("*", ".", "-")]
  # Kinds of aa or nt
  n_elements <- length(unique(seq))
  if (n_elements == 0) {
    var = 0               # Set variability of gap site as 0
    # var = NA            # Set variability of gap site as NA
  } else {
    # Frequency of the top aa or nt
    freq = data.frame(table(seq))
    top_freq = max(freq$Freq / sum(freq$Freq))
    # Calculate variability 
    var = ifelse(n_elements > 0, n_elements / top_freq, NA)
  }
  return(var)
}

## Function: Calculate variability indexes for all sites
calculate_variabilities <- function(seqs) {
  n_sites <- nchar(seqs[1])
  sites <- lapply(1:n_sites, function(i) substr(seqs, i, i))
  variabilities <- sapply(sites, calculate_variability)
  return(variabilities)
}


# Function: Calculate variability indexes in all sites for all samples
calculate_vdj_variabilities <- function(vdj_list) {
  lapply(vdj_list, function(vdj) {
    sample_name <- unique(vdj$sample)
    v_index  <- data.frame(variability = calculate_variabilities(vdj$imgt_aa)) %>%
      rownames_to_column("position") %>%
      mutate(position = as.numeric(position)) %>%
      mutate(region = case_when(position < 27 ~ "FR1",
                                position < 39 ~ "CDR1",
                                position < 56 ~ "FR2",
                                position < 66 ~ "CDR2",
                                position < 105 ~ "FR3")) %>%
      mutate(sample = sample_name)
  }) %>% bind_rows()
}


v_batrian <- v_gdna %>% select(type, imgt_aa) %>% mutate(species = "Camelus bactrianus") 
vdj <- v_camelids %>% subset(molecule=="gDNA") %>% 
  bind_rows(v_batrian) %>%
  select(species, type, imgt_aa)

vh_vdj_list <- vdj %>% subset(type == "VH") %>% group_split(species)
vhh_vdj_list <- vdj %>% subset(type == "VHH") %>% group_split(species)
vh_vdj_variabilities <- calculate_vdj_variabilities(vh_vdj_list) %>% mutate(type = "VH")
vhh_vdj_variabilities <- calculate_vdj_variabilities(vhh_vdj_list) %>% mutate(type = "VHH")
vdj_variabilities <- bind_rows(vhh_vdj_variabilities, vh_vdj_variabilities) 



  
p_var <- vdj_variabilities %>%
  mutate(variability = ifelse(type == "VH", -variability, variability)) %>%
  ggplot() +
  geom_hline(yintercept=0, color = "black") +
  geom_bar(mapping = aes(x = position, y = variability, fill = region),
           stat = "identity", width = 0.6) +
  scale_x_continuous(breaks = seq(0, 105, by = 1)) +
  facet_grid(species ~ .) +
  theme_classic() +

  # annotate(geom = "text", x = 5, y = -10, label = "VH") +
  # annotate(geom = "text", x = 5, y = 10, label = "VHH") +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, colour = "black", hjust=0.95,vjust=0.3),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Position", y = "Variability index", fill = "Region") 


pdf("results/camelids/variability_plot.pdf", width = 8, height = 5)
p_var
dev.off()






p_diff <- vdj_variabilities %>%
  mutate(variability = ifelse(type == "VH", -variability, variability)) %>%
  group_by(position, region, species) %>%
  summarise(diff = sum(variability)) %>%
  mutate(change = ifelse(diff > 0, "up", "down")) %>%
  mutate(diff = abs(diff)) %>%
  ggplot() +
  geom_tile(aes(x = position, y = "1", fill = change, alpha = diff), width = 0.9) +
  # geom_text(aes(x = position, y = "1", label = round(diff,1))) +
  facet_grid(species ~ . ) +
  scale_fill_manual(values = c("#34bf49", "#ff4c4c", "#caccd1")) +
  scale_x_continuous(breaks = seq(0, 105, by = 1))  +
  facet_grid(species ~ . ) +
  theme_classic() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle=90, colour="black", hjust=0.95, vjust=0.3),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(caption = paste0("chisq, p_value=", chisq_test$p.value))
p_diff

pdf("results/camelids/variability_plot_diff.pdf", width = 8, height = 2)
p_diff
dev.off()



changes <- vdj_variabilities %>%
  mutate(variability = ifelse(type == "VH", -variability, variability)) %>%
  group_by(position, region, species) %>%
  summarise(diff = sum(variability)) %>%
  mutate(change = ifelse(diff > 0, "up", "down")) %>%
  group_by(species, change) %>% summarise(count = n())

chisq_test <- vdj_variabilities %>%
  mutate(variability = ifelse(type == "VH", -variability, variability)) %>%
  group_by(position, region, species) %>%
  summarise(diff = sum(variability)) %>%
  mutate(change = ifelse(diff > 0, "up", "down")) %>%
  group_by(species, change) %>% summarise(count = n()) %>%
  spread(change, count) %>% 
  column_to_rownames("species") %>%
  chisq.test()
chisq_test

p_changes <- changes %>% 
  group_by(species) %>% mutate(prop = count / sum(count)) %>%
  ggplot(aes(x=species, y = prop, fill = change)) +
  scale_fill_manual(values = c("#34bf49", "#ff4c4c", "#caccd1")) +
  geom_bar(stat = "identity") + 
  theme_classic() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle=90, colour="black", hjust=0.95, vjust=0.3),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(caption = paste("chisq_test,", chisq_test$p.value))
p_changes

pdf("results/camelids/variability_plot_changes.pdf", width = 3, height = 6)
p_changes
dev.off()
