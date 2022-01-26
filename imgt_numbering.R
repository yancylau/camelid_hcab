#
# Fucntions for IMGT unique numbering



#### IMGT numbering for FRs
# Note: seq is a vector with FR sequences
# 内部实质都是向量化运算
fr_numbering_aa <- function(seq, region) {
  seq_len = nchar(seq)

  # MGT numbering for FR1
  if (region == "FR1") {
    imgt_aa = case_when(
      seq_len == 25 ~ paste0(substr(seq, 1, 9), ".", substr(seq, 10, 25)),
      seq_len == 26 ~ seq,
      TRUE ~ ""
    )
    return(imgt_aa)
  }
  # MGT numbering for FR2
  if (region == "FR2") {
    imgt_aa = ifelse(seq_len == 17, seq, "")
    return(imgt_aa)
  }
  # MGT numbering for FR3
  if (region == "FR3") {
    imgt_aa = case_when(
      seq_len == 38 ~ paste0(substr(seq, 1, 7), ".", substr(seq, 8, 38)),
      seq_len == 39 ~ seq,
      TRUE ~ ""
    )
    return(imgt_aa)
  }
}


fr_numbering_nt <- function(seq, region) {
  seq_len = nchar(seq)
  
  # MGT numbering for FR1
  if (region == "FR1") {
    imgt_nt = case_when(
      seq_len == 75 ~ paste0(substr(seq, 1, 27), "...", substr(seq, 28, 75)),
      seq_len == 78 ~ seq,
      TRUE ~ ""
    )
    return(imgt_nt)
  }
  # MGT numbering for FR2
  if (region == "FR2") {
    imgt_nt = ifelse(seq_len == 51, seq, "")
    return(imgt_nt)
  }
  # MGT numbering for FR3
  if (region == "FR3") {
    imgt_nt = case_when(
      seq_len == 114 ~ paste0(substr(seq, 1, 21), "...", substr(seq, 22, 114)),
      seq_len == 117 ~ seq,
      TRUE ~ ""
    )
    return(imgt_nt)
  }
}


#### IMGT numbering for CDR3
# Note: Seqs is vector with CDR sequences
cdr_numbering_aa <- function(seqs, region) {
  if (region == "CDR1") {
    imgt_length = 12
    gaps = c(7, 6, 8, 5, 9, 4, 10)
  } 
  if (region == "CDR2") {
    imgt_length = 10
    gaps = c(6, 5, 7, 4, 8, 3, 9, 2, 10, 1)
  } 
  if (region == "CDR3") {
    imgt_length = 13
    gaps = c(7, 8, 6, 9, 5, 10, 4, 11)
  }
  
  lapply(seqs, function(seq) {
    seq_len = nchar(seq)
    
    n_gaps = imgt_length - seq_len
    gap_sites = gaps[1:n_gaps]
    gaps_start = min(gap_sites)
    # Insert gaps to seq
    paste0(substr(seq, 1, gaps_start - 1), paste0(rep(".", n_gaps), collapse = ""), substr(seq, gaps_start , seq_len), collapse = "")
  }) %>% unlist()
}



cdr_numbering_nt <- function(seqs, region) {
  if (region == "CDR1") {
    imgt_length = 12 
    gaps = c(7, 6, 8, 5, 9, 4, 10)
  } 
  if (region == "CDR2") {
    imgt_length = 10
    gaps = c(6, 5, 7, 4, 8, 3, 9, 2, 10, 1)
  } 
  if (region == "CDR3") {
    imgt_length = 13
    gaps = c(7, 8, 6, 9, 5, 10, 4, 11)
  }
  
  lapply(seqs, function(seq) {
    seq_len = nchar(seq)
    
    n_gaps = imgt_length - seq_len / 3
    gap_sites = gaps[1:n_gaps]
    gaps_start = (min(gap_sites) - 1) * 3 + 1
    # Insert gaps to seq
    paste0(substr(seq, 1, gaps_start - 1), paste0(rep(".", n_gaps * 3), collapse = ""), substr(seq, gaps_start, seq_len), collapse = "")
  }) %>% unlist()
}


# cdr_numbering <- function(seq, region) {
#   if (region == "CDR1") {
#     imgt_length = 12
#     gaps = c(7, 6, 8, 5, 9, 4, 10)
#   } 
#   if (region == "CDR2") {
#     imgt_length = 10
#     gaps = c(6, 5, 7, 4, 8, 3, 9, 2, 10, 1)
#   } 
#   if (region == "CDR3") {
#     imgt_length = 13
#     gaps = c(7, 8, 6, 9, 5, 10, 4, 11)
#   }
#   
#   seq_len = nchar(seq)
#   
#   n_gaps = imgt_length - seq_len
#   gap_sites = gaps[1:n_gaps]
#   gaps_start = min(gap_sites)
#   # Insert gaps to seq
#   paste0(substr(seq, 1, gaps_start), paste0(rep(".", n_gaps), collapse = ""), substr(seq, gaps_start, seq_len), collapse = "")
# }







