# Function to calculate variability in a single site
cal_var <- function(seq) {
  # Remove gaps or in this site
  seq <- seq[! seq %in% c("*", ".", "-")]
  # Kinds of aa or nt
  n_elements <- length(unique(seq))
  if (n_elements == 0) {
    var = 0        # Set variability of gap site as 0
    # var = NA     # Set variability of gap site as NA
  } else {
    # Frequency of the top aa or nt
    freq = data.frame(table(seq))
    top_freq = max(freq$Freq / sum(freq$Freq))
    # Calculate variability 
    var = ifelse(n_elements > 0, n_elements / top_freq, NA)
  }
  return(var)
}
# Function to calculate variability of all sites
cal_var_index <- function(seqs) {
  n_sites <- nchar(seqs[1])
  sites <- lapply(1:n_sites, function(i) substr(seqs, i, i))
  variabilities <- sapply(sites, cal_var)
  return(variabilities)
}
