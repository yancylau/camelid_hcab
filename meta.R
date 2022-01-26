


codons <- data.frame(codon = c("GCC", "GCG", "GCU", "GCA", "AGA", "CGG", "AGG", "CGA", 
                               "CGC", "CGU", "AAC", "AAU", "GAC", "GAU", "UGC", "UGU", 
                               "CAA", "CAG", "GAG", "GAA", "GGC", "GGU", "GGA", "GGG", 
                               "CAC", "CAU", "AUA", "AUC", "AUU", "CUG", "CUA", "UUA", 
                               "CUU", "UUG", "CUC", "AAA", "AAG", "AUG", "UUC", "UUU", 
                               "CCG", "CCC", "CCU", "CCA", "AGC", "UCG", "UCU", "UCA", 
                               "UCC", "AGU", "UAG", "UAA", "UGA", "ACA", "ACC", "ACG", 
                               "ACU", "UGG", "UAU", "UAC", "GUA", "GUG", "GUU", "GUC"),
                    aa = c("A", "A", "A", "A", "R", "R", "R", "R", 
                           "R", "R", "N", "N", "D", "D", "C", "C", 
                           "Q", "Q", "E", "E", "G", "G", "G", "G", 
                           "H", "H", "I", "I", "I", "L", "L", "L", 
                           "L", "L", "L", "K", "K", "M", "F", "F",
                           "P", "P", "P", "P", "S", "S", "S", "S", 
                           "S", "S", "*", "*", "*", "T", "T", "T", 
                           "T", "W", "Y", "Y", "V", "V", "V", "V")) %>% 
  mutate(codon = str_replace_all(codon, "U", "T"))
