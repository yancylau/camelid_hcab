# Function to calculate allele frequence of VHH and VH
# Function to compare allele frequence of VHH and VH

# Argument: dataframe with "type" and "frs" columns
#   type: VH/VHH
#   frs: imgt numbered FRs (FR1+FR2+FR3)
# Return: dataframe




source("R/chemistry.R")

# Color pallate
aas <- c("A", "C", "D", "E", "F", 
         "G", "H", "I", "K", "L",
         "M", "N", "P", "Q", "R", 
         "S", "T", "V", "W", "Y")


# colors <- c("#6E6157", "#BCBAAD", "#A5A63A", "#ABC340", "#12575F",
#               "#84A49A", "#657790", "#8FADCC", "#9DCECE", "#258EC6",
#               "#004570", "#3F8AA3", "#DD6338", "#ECC968", "#CE8734",
#               "#DF8E92", "#BC2374", "#7B5C9D", "#4B2974", "#C82424",
#               "#76191C", "#060001", "#7A7F83")

colors <- c("#2bb594","#15cddc","#31a7b7","#4276b3","#714ca3",
            "#b74fa2","#f04f99","#ff4b48","#ff733f","#ff9936",
            "#037171","#03312e","#0e0f19","#89608e","#e9d6ec",
            "#b2d3a8","#ede5a6","#36f1cd","#3bc14a","#c47ac0")


# Allele frequency --------------------------------------------------------
site_dict <- tibble(pos = c(1:82), site = c(1:26, 39:55, 66:104))

calc_allele_freq <- function(frs) {
  
  # Allete freq of all sites
  allele_freq <- function(frs) {
    seq_len = nchar(frs$frs)[1]
    c(1:seq_len) %>% map(function(pos) {
      frs %>% 
        mutate(composite = substr(frs, pos, pos)) %>% 
        group_by(composite) %>% summarise(count = n()) %>%
        mutate(prop = count / sum(count))%>%
        mutate(pos = pos)
    }) %>% bind_rows()
  }
  
  prop_vh <- frs %>% subset(type == "VH") %>% allele_freq() %>% mutate(type = "VH")
  prop_vhh <- frs %>% subset(type == "VHH") %>% allele_freq() %>% mutate(type = "VHH")
  bind_rows(prop_vh, prop_vhh) %>%
    left_join(site_dict, by = "pos") %>%
    select(-pos)
}



# Find differential sites ----------------------------------------------------
chisq_test <- function(frs) {
  seq_len = nchar(frs$frs)[1]
  c(1:seq_len) %>% map(function(pos) {
    frs %>% mutate(composite = substr(frs, pos, pos)) %>%
      select(type, composite) %>% 
      table() %>% chisq.test() %>% broom::tidy() %>% 
      set_names(str_replace_all(names(.), "\\.", "_")) %>%
      mutate(pos = pos)
  }) %>% bind_rows() %>%
    left_join(site_dict, by = "pos") %>% select(-pos)
} 



get_diff_sites <- function(frs) {
  aa_freqs <- frs %>% calc_allele_freq() 
  
  # Chi-square test
  chisq <- frs %>% chisq_test() 
  
  # Subset diff sites
  diff_sites <- chisq %>% subset(p_value < 0.05) %>% pull(site)
  # Subset major site
  major_sites <- aa_freqs %>% 
    subset(prop > 0.8) %>% 
    arrange(site) %>% 
    group_by(site, composite) %>% 
    summarise(count = n()) %>%
    subset(count < 2) %>% pull(site) %>% unique()
  
  intersect(major_sites, diff_sites)
}



chemical_chisq_test <- function(frs) {
  seq_len = nchar(frs$frs)[1]
  
  lapply(1:seq_len, function(pos){
    chems <- frs %>% 
      mutate(composite = substr(frs, pos, pos)) %>%
      left_join(hydropathy, by = "composite") %>%
      left_join(volume, by = "composite") %>%
      left_join(chemical, by = "composite") %>%
      left_join(charge, by = "composite") %>% 
      left_join(hydrogen, by = "composite") %>% 
      left_join(polarity, by = "composite") %>%
      select(type, hydropathy, volume, chemical, charge, hydrogen, polarity) %>% 
      gather(key = "chemical", value = "property", -type) 
      
    
    chemicals <- c("hydropathy", "volume", "chemical", "charge", "hydrogen", "polarity")
    lapply(chemicals, function(chem) {
      n_levels <- chems %>% subset(chemical == chem) %>% select(property) %>% n_distinct()
      if(n_levels < 2) {return(NULL)}
      
      chems %>% subset(chemical == chem) %>%
        select(type, property) %>% table() %>% chisq.test() %>% broom::tidy() %>%
        mutate(chemical = chem)
    }) %>% bind_rows() %>% mutate(pos = pos) 
  }) %>% bind_rows() %>%
    left_join(site_dict, by = "pos") %>% select(-pos)
}



# Plots -------------------------------------------------------------------
plot_composite <- function(aa_freqs, sites = c(1:26, 39:55, 66:104)) {
  df <- aa_freqs %>%
    subset(site %in% sites) %>%
    mutate(region = case_when(
      site < 27 ~ "FR1",
      site < 56 ~ "FR2",
      site < 105 ~ "FR3"
    )) %>%
    mutate(composite = str_replace(composite, "\\.", "z")) %>%
    mutate(site = factor(site, levels = c(1:26, 39:55, 66:104))) %>%
    mutate(type = factor(type, levels = c("VHH", "VH"))) %>%
    mutate(composite = fct_reorder(composite, prop, .fun=sum))
  
  ggplot(df, aes(x = site, y = prop, group = count, fill = composite)) +
    geom_bar(stat = "identity", position = "fill", width = 0.8) +
    theme_classic() +
    theme(strip.background = element_blank(),
          # axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2, colour = "black"),
          axis.text = element_text(color = "black"),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.position = "bottom") +
    facet_grid(type ~ region, scales = "free", space = "free") +
    labs(x = "", y = "Proportion", fill = "Amino acid") +
    guides(fill=guide_legend(nrow = 1, byrow = TRUE)) +
    scale_fill_manual(values = colors, labels = aas)
}


plot_chemical <- function(aa_freqs, sites = c(1:26, 39:55, 66:104)) {
  df <- aa_freqs %>%
    subset(site %in% sites) %>%
    left_join(hydropathy, by = "composite") %>%
    left_join(volume, by = "composite") %>%
    left_join(chemical, by = "composite") %>%
    left_join(charge, by = "composite") %>% 
    left_join(hydrogen, by = "composite") %>% 
    left_join(polarity, by = "composite") %>%
    gather("chemical", "property", -c(composite,count, prop, type, site)) %>%
    group_by(site, chemical, type, property) %>% summarise(count = sum(count)) %>%
    group_by(site, chemical, type) %>% mutate(prop = count/sum(count)) %>%
    mutate(region = case_when(
      site < 27 ~ "FR1",
      site < 56 ~ "FR2",
      site < 105 ~ "FR3"
    )) %>%
    mutate(site = factor(site, levels = c(1:26, 39:55, 66:104))) %>%
    mutate(type = factor(type, levels = c("VHH", "VH"))) 
  
  df_list <- df %>% ungroup() %>% group_split(chemical)
  df_list %>% lapply(function(df){
    df %>% 
      mutate(property = fct_reorder(property, prop, .fun=sum)) %>%
      ggplot(aes(x = site, y = prop, fill = property)) +
      geom_bar(stat = "identity", position = "fill", width = 0.8) +
      theme_classic() +
      theme(strip.background = element_blank(),
            # axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2, colour = "black"),
            axis.text = element_text(color = "black"),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(colour = "black"),
            legend.position = "bottom") +
      facet_grid(type ~ chemical + region, scales = "free", space = "free") +
      labs(x = "", y = "Proportion") +
      guides(fill=guide_legend(nrow = 1, byrow = TRUE))
  })
}


# Main function -----------------------------------------------------------
composite_test <- function(frs, path) {
  chisq_sites <- frs %>% chisq_test() 
  write_csv(chisq_sites, path)
}

composite_plot <- function(frs, path, ...) {
  aa_freqs <- frs %>% calc_allele_freq() 
  diff_sites <- get_diff_sites(frs)
  p_diff_sites <- plot_composite(aa_freqs, diff_sites)
  
  pdf(path, ...)
  print(p_diff_sites)
  dev.off()
}


chemicals_test <- function(frs, path) {
  chisq_chemicals <- frs %>% chemical_chisq_test() 
  write_csv(chisq_chemicals, path)
}

chemicals_plot <- function(frs, path, ...) {
  aa_freqs <- frs %>% calc_allele_freq() 
  diff_sites <- get_diff_sites(frs)
  p_diff_chems <- plot_chemical(aa_freqs, diff_sites)
  
  pdf(path, ...)
  print(p_diff_chems)
  dev.off()
}





