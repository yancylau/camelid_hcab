# Calculate variability index on all sites




# Variability index on a certain site
calculate_variability <- function(aas, meth = "classical") {
  aas <- aas[! aas %in% c("*", ".", "-")]   # Remove gaps
  n_elements <- length(unique(aas))
  
  if (n_elements == 0) {return(0) }    # Set gap as 0 or NA
  
  freq <- data.frame(table(aas))     
  if (meth == "classiical") {
    top_freq <- max(freq$Freq / sum(freq$Freq))
    v_index <- n_elements / top_freq
  } else {
    top_freq <- max(freq$Freq / sum(freq$Freq))
    v_index <- top_freq
  }
  
  return(v_index)
}



calculate_variability <- function(seq) {   # Remove gaps
  seq <- seq[! seq %in% c("*", ".", "-")]
  n_elements <- length(unique(seq))
  
  if (n_elements == 0) {
    return(0)                              # Set gap as 0 or NA
  } else {
    freq = data.frame(table(seq))
    top_freq = max(freq$Freq / sum(freq$Freq))
    # Calculate variability 
    var = ifelse(n_elements > 0, n_elements / top_freq, NA)
  }
  return(var)
}



# Calculate variability indexes of all sites
calculate_variabilities <- function(seqs) {
  n_sites <- nchar(seqs[1])
  sites <- lapply(1:n_sites, function(i) substr(seqs, i, i))
  variabilities <- map_dbl(sites, calculate_variability)
  return(variabilities)
}




# Make variability index table
make_indexes_table <- function(df) {
  vhh_seqs <- df %>% subset(type == "VHH") %>% select(imgt_aa) %>% pull() 
  vh_seqs <- df %>% subset(type == "VH") %>% select(imgt_aa) %>% pull() 

  vhh_v_indexes <- calculate_variabilities(vhh_seqs) %>% 
    data.frame() %>% set_names("variability") %>% 
    rownames_to_column("position") %>%
    mutate(type = "VHH")
  vh_v_indexes <- calculate_variabilities(vh_seqs) %>% 
    data.frame() %>% set_names("variability") %>% 
    rownames_to_column("position") %>%
    mutate(type = "VH")

  bind_rows(vhh_v_indexes, vh_v_indexes) %>%
    mutate(position = as.numeric(position)) %>%
    # mutate(variability = ifelse(type == "VH", -variability, variability)) %>%
    mutate(region = case_when(position < 27 ~ "FR1",
                              position < 39 ~ "CDR1",
                              position < 56 ~ "FR2",
                              position < 66 ~ "CDR2",
                              position < 105 ~ "FR3")) 
}



# Function to plot variability index
plot_variability <- function(v_indexes) {
  v_indexes %>% 
    mutate(variability = ifelse(type == "VH", -variability, variability)) %>%
    mutate(type = factor(type, levels = c("VHH", "VH"))) %>%
    mutate(position = factor(position)) %>%
    ggplot(aes(x = position, y = variability, fill = position, color = position)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_x_discrete(breaks = seq(1, 104, by = 1)) +
    theme_classic() +
    # geom_hline(yintercept=0, color = "black") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.3)) +
    labs(x = "Position", y = "Variability index", fill = "Region") +
    facet_grid(type ~.) 
}


plot_variability_diff <- function(v_indexes) {
  v_indexes %>%
    spread(type, variability) %>%
    mutate(diff = VHH - VH) %>%
    mutate(change = case_when(diff > 0 ~ "Up",
                              diff < 0 ~ "Down",
                              diff == 0 ~ "Unchange")) %>%
    mutate(change = factor(change, levels = c("Down", "Up", "Unchange"))) %>%
    mutate(value = ifelse(change == "Unchange", 1, abs(diff))) %>%
    ggplot() +
    geom_tile(aes(x = position, y = "1", fill = change, alpha = value), width = 0.7) +
    scale_fill_manual(values = c("#34bf49", "#ff4c4c", "grey")) +
    scale_x_continuous(breaks = seq(0, 105, by = 1))  +
    theme_classic() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 90, colour = "black", hjust = 0.95, vjust = 0.3),
          axis.text.y = element_text(colour = "black"),
          axis.ticks = element_line(color = "black"))
}


plot_variability_change <- function(v_indexes) {
  v_indexes %>%
    spread(type, variability) %>%
    mutate(diff = VHH - VH) %>%
    mutate(change = case_when(diff > 0 ~ "Up",
                              diff < 0 ~ "Down",
                              diff == 0 ~ "Unchange")) %>%
    mutate(change = factor(change, levels = c("Down", "Up", "Unchange"))) %>%
    mutate(value = abs(diff)) %>%
    ggplot(aes(x = position, y = value, fill = change)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = c("#34bf49", "#ff4c4c", "grey")) +
    scale_x_continuous(breaks = seq(0, 105, by = 1))  +
    theme_classic() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 90, colour = "black", hjust = 0.95, vjust = 0.3),
          axis.text.y = element_text(colour = "black"),
          axis.ticks = element_line(color = "black"))
}



to_pdf <- function(expr, filename, ..., verbose=TRUE) {
  if (verbose)
    cat(sprintf("Creating %s\n", filename))
  pdf(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

# to_pdf(plot_variability(v_indexes), "results/test2.pdf", width=6, height=4)



