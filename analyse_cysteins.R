#

colors <- c("#2bb594","#15cddc","#31a7b7","#4276b3","#714ca3",
            "#b74fa2","#f04f99","#ff4b48","#ff733f","#ff9936",
            "#037171","#03312e","#0e0f19","#89608e","#e9d6ec",
            "#b2d3a8","#ede5a6","#36f1cd","#3bc14a","#c47ac0")




# Summary cysteins --------------------------------------------------------
summary_cysteins <- function(v) {
  v %>% 
    # rename(full = imgt_aa) %>%
    # gather(key = "region", value = "seq", 
    #        -c(id, family, species, molecule, subgene, type, from)) %>%
    mutate(n_cys = str_count(imgt_aa, "C")) %>%
    select(-id, -imgt_aa) %>% group_by_all() %>% summarise(count = n()) %>%
    group_by(family, species, subgene, type, from) %>%
    mutate(total = sum(count)) %>%
    mutate(prop = count / total) 
}




# Find cysteins -----------------------------------------------------------
get_c_positions <- function(v) {
  position_list <- v %>% pull(imgt_aa) %>% str_locate_all("C")
  
  lapply(seq_along(position_list), function(i) {
    id = v$id[i]
    pos = position_list[[i]][,1]
    data.frame(id = id, position = pos)
  }) %>% bind_rows() 
}





cystein_pattern <- function(v) {
  get_c_positions(v) %>% 
    mutate(region = case_when(
      position < 27 ~ "FR1",
      position < 39 ~ "CDR1",
      position < 56 ~ "FR2",
      position < 66 ~ "CDR2",
      position < 105 ~ "FR3"
    )) %>%
    group_by(id) %>% add_count(name = "n_cys") %>%
    mutate(cystein = str_c(region, "(", position, ")")) %>% select(-position, -region) %>%
    mutate(count = 1) %>% spread(cystein, count, fill = 0) %>%
    ungroup() %>% select(-id) %>%
    group_by_all() %>% summarise(count = n()) %>%
    ungroup() %>% mutate(prop = count / sum(count))
}



summary_c_positions <- function(v) {
  position_list <- v %>% pull(imgt_aa) %>% str_locate_all("C")
  
  cys_positions <- lapply(seq_along(position_list), function(i) {
    id = v$id[i]
    pos = position_list[[i]][,1]
    data.frame(id = id, position = pos)
  }) %>% bind_rows() 
  
  v %>% 
    select(id, family, species, molecule, subgene, type) %>%
    group_by(family, species, molecule, subgene, type) %>% 
    add_count(name = "total") %>%
    left_join(cys_positions, by = "id") %>%
    select(-id) %>% group_by_all() %>% summarise(count = n()) %>%
    mutate(prop = count/total)
}




# Compare cysteins --------------------------------------------------------
# Chi square test on cystein groups

test_cysteins <- function(n_cys) {
  n_cys %>% 
    subset(total > 5) %>%
    mutate(group = ifelse(n_cys > 2, ">=3", "2")) %>%
    group_by(family, species, type, group, total) %>% 
    summarise(count = sum(count)) %>%
    ungroup() %>% select(species, type, group, count) %>%
    complete(species, type, nesting(group), fill = list(count = 0)) %>%
    spread(group, count) %>% select(-type) %>% 
    nest_by(species) %>% mutate(model = list(chisq.test(data))) %>%
    summarise(tidy(model))
}


# Plot cystein comparison -------------------------------------------------
plot_cystein <- function(n_cys) {
  df <- n_cys %>% subset(total > 5) %>%
    mutate(group = ifelse(n_cys > 2, ">=3", "2")) %>%
    group_by(family, species, subgene, type, from, group, total) %>% summarise(count = sum(count)) %>%
    spread(group, count, fill = 0) %>% gather("group", "count", -c(family, species, subgene, type, from, total)) %>%
    mutate(prop = round(count/total, 4))  %>%
    mutate(type = factor(type, levels = c("VHH", "VH"))) %>%
    mutate(species = reorder_within(species, prop, list(family, subgene, type)))

  ggplot(df, aes(x = species, y = prop)) +
  # ggplot(df, aes(x = reorder_within(species, prop, list(family, subgene, type)), y = prop)) +
    geom_bar(stat = "identity", aes(fill = group)) +
    scale_x_reordered() +
    facet_grid(. ~ family + subgene + type, scales = "free", space = "free") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, 
                                     colour = "black", 
                                     hjust = 0.95, 
                                     vjust = 0.2)) +
    labs(x = "", y = "Proportion", fill = "Species") 
}




# Plot cystein distribution -----------------------------------------------
plot_cystein_distribution <- function(c_positions) {
    df <- c_positions %>%
    subset(total > 5) %>%
    mutate(region = case_when(
      position < 27 ~ "FR1",
      position < 39 ~ "CDR1",
      position < 56 ~ "FR2",
      position < 66 ~ "CDR2",
      position < 105 ~ "FR3"
    )) %>%
    mutate(positon = factor(position, levels = c(1:104))) %>%
    mutate(group = ifelse(position %in% c(23, 104), "Canonical", as.character(position))) %>%
    # subset(family == "camelids") %>%
    subset(subgene == "IGHV3") %>%
    mutate(type = factor(type, levels = c("VHH", "VH")))
  
  # IGHV3 + non-IGHV3
  ggplot(df, aes(x = position, y = 0, fill = group)) +
    geom_rect(mapping=aes(xmin=0, xmax=26.5, ymin=-1, ymax=1), fill=NA, color = "black") + 
    geom_rect(mapping=aes(xmin=26.5, xmax=38.5, ymin=-1, ymax=1), fill=NA, color = "black") +
    geom_rect(mapping=aes(xmin=38.5, xmax=55.5, ymin=-1, ymax=1), fill=NA, color = "black") +
    geom_rect(mapping=aes(xmin=55.5, xmax=66.5, ymin=-1, ymax=1), fill=NA, color = "black") +
    geom_rect(mapping=aes(xmin=66.5, xmax=105, ymin=-1, ymax=1), fill=NA, color = "black") +
    geom_point(aes(size = prop), shape=21) +
    scale_size(range = c(.1, 10)) +
    theme_classic() + 
    scale_fill_manual(values = colors) +
    theme(
      axis.line.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()) +
    scale_x_continuous(breaks = c(1:105), labels = c(1:105))+
    facet_grid(~ family + type + species ~ .) +
    theme(axis.text.x = element_text(angle = 90, colour = "black", hjust = 0.95, vjust = 0.2))
}









