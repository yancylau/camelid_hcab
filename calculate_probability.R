#!/usr/bin/env Rscript
#
#
#
suppressMessages(library(tidyverse))
load("data/merge/vdj.Rda")

seqs = vdj$imgt_aa

aas <- c(
  'A', 'R', 'N', 'D', 'C', 
  'Q', 'E', 'G', 'H', 'I', 
  'L', 'K', 'M', 'F', 'P', 
  'S', 'T', 'W', 'Y', 'V'
)
get_letter_matrix <- function(seqs){
  # Ensure kmers are the same length characters 
  seq_len = sapply(seqs, nchar)
  num_pos = seq_len[1]
  if(! all(seq_len == num_pos)) stop('Sequences length are not identical')
  # Construct matrix of letters
  split = unlist( sapply(seqs, function(seq){strsplit(seq, '')}) )
  t( matrix(split, seq_len, length(split)/num_pos) )
}
aa_mat <- get_letter_matrix(seqs)
namespace <- intersect(as.character(aa_mat), aas)

# Construct PWM
ratio_mat = apply(
  aa_mat, 2, function(position){
    # Get frequencies 
    t = table(position)
    # Match to aa
    ind = match(namespace, names(t))
    # Create column
    col = t[ind]
    col[is.na(col)] = 0
    names(col) = namespace
    # Do relative frequencies
    col = col / sum(col)
    col
  }
)


 
# Load data
load("data/clustering/atleast_2/vdj.Rda")
seqs <- vdj$cdr1
seqs <- vdj$imgt_aa



.AA_NAMESPACE = function() c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
letterMatrix <- function(input){
  # Ensure kmers are the same length characters 
  seq.len = sapply(input, nchar)
  num_pos = seq.len[1]
  if(! all(seq.len == num_pos)) stop('Sequences in alignment must have identical lengths')
  
  # Construct matrix of letters
  split = unlist( sapply(input, function(seq){strsplit(seq, '')}) )
  
  t( matrix(split, seq.len, length(split)/num_pos) )
}
letter_mat <- letterMatrix(seqs)
namespace <- intersect(as.character(letter_mat), .AA_NAMESPACE())

# Construct PWM
pfm_mat = apply(letter_mat, 2, function(pos.data){
  # Get frequencies 
  t = table(pos.data)
  # Match to aa
  ind = match(namespace, names(t))
  # Create column
  col = t[ind]
  col[is.na(col)] = 0
  names(col) = namespace
  # Do relative frequencies
  col = col / sum(col)
  col
})




