



# Detect gaps and deletes
detect_gaps <- function(v_gdna, v_cdna) {
  gaps_gdna <- str_locate_all(v_gdna, "[.]") %>% unlist() %>% unique() 
  gaps_cdna <- str_locate_all(v_cdna, "[.]") %>% unlist() %>% unique()     
  # Inserts
  inserts = setdiff(gaps_gdna, gaps_cdna)   # Gaps only in germline v
  if (length(inserts) > 0) {
    df_inserts = data.frame(site = inserts, variant = "insert")
  } else{
    df_inserts = NULL
  }
  
  # Deletes
  deletes = setdiff(gaps_cdna, gaps_gdna)   # Gaps only  in rearranged v
  if (length(deletes) > 0) {
    df_deletes = data.frame(site = deletes, variant = "delete")
  } else{
    df_deletes = NULL
  }
  # df_inserts = ifelse(length(inserts) > 0, data.frame(site = inserts, variant = "insert"), NULL)
  # df_deletes = ifelse(length(deletes) > 0, data.frame(site = deletes, variant = "insert"), NULL)
  
  df_gaps <- bind_rows(df_inserts, df_deletes)
}



# Delete mismatches
detect_mismatches <- function(v_gdna, v_cdna) {
  v = vector(mode = "numeric", length = 0)
  for (i in c(1:104)) {
    if(str_sub(v_gdna, i, i) != str_sub(v_cdna, i, i)) v = c(v, i)
  }
  
  if (length(v) > 0){
    df_mismatches <- data.frame(site = v, variant = "mismatches")
  } else {
    df_mismatches <- NULL
  }
  
  return(df_mismatches)
}





# Plot mismatches
