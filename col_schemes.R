
# Adapoted from ggseqlog



# List color schemes available in ggseqlogo
list_col_schemes <- function(v=T){
  
  col_schemes = c('auto', 'chemistry', 'chemistry2','hydrophobicity', 'nucleotide', 'nucleotide2',
                  'base_pairing', 'clustalx', 'taylor')
  if(!v) return(col_schemes)
  message('Available ggseqlogo color schemes:')
  for(f in col_schemes) message('\t', f)
}



# Get color scheme
get_col_scheme = function(col_scheme, seq_type='auto'){
  
  # Check if user-defined color scheme
  if(is.data.frame(col_scheme)){
    if(!'ggseqlogo_cs' %in% class(col_scheme)) 
      stop('Colour scheme must be generated using "make_col_scheme" function')
    return(col_scheme)
  }
  
  # Get ambigious colour scheme
  # col_scheme = match.arg(col_scheme, list_col_schemes(F))
  
  # Get default color scheme for sequence type
  if(col_scheme == 'auto'){
    if(seq_type == 'auto') stop('"col_scheme" and "seq_type" cannot both be "auto"')
    
    col_scheme = switch(tolower(seq_type), aa = 'chemistry', 
                        dna = 'nucleotide', rna = 'nucleotide', 
                        other='nucleotide')
    
  }
  
  
  # Pick from default color schemes
  cs = switch(col_scheme, 
              
              # 'Hydropathy'
              hydropathy = data.frame(
                letter = c(
                  c("A", "C", "I", "L", "M", "F", "W", "V"), 
                  c("G", "H", "P", "S", "T", "Y"), 
                  c("R", "N", "D", "Q", "E", "K")
                ),
                group = c(
                  rep("Hydrophobic", 8), 
                  rep("Neutral", 6), 
                  rep("Hydrophilic", 6)
                ),
                col = c(
                  rep("0099e5", 8), 
                  rep("ff4c4c", 6), 
                  rep("34bf49", 6)
                ),
                stringsAsFactors = F
              ),
              
              
              # 'Chemical'
              chemical = data.frame(
                letter = c(
                  c("A", "G", "I", "L", "P", "V"),
                  c("F", "W", "Y"),
                  c("C", "M"),
                  c("S", "T"),
                  c("R", "H", "K"),
                  c("D", "E"),
                  c("N", "Q")
                ),
                group = c(
                  rep("aliphatic ", 6), 
                  rep("aromatic", 3), 
                  rep("sulfur", 2),
                  rep("hydroxyl", 2), 
                  rep("basic", 3), 
                  rep("acidic", 2),
                  rep("amide", 2)
                ),
                col = c(
                  rep("#222c37", 6), 
                  rep("#00cccc", 3), 
                  rep("#fff600", 2),
                  rep("#ff0066", 2), 
                  rep("#19e3b1", 3), 
                  rep("#ff7f33", 2),
                  rep("#b83c82", 2)
                ),
                stringsAsFactors = F
              ),

              # Hydrophobicity index (PMID: 7108955) from -4.5 to 4.5
              hydrophobicity = data.frame(
                letter = c('I', 'V', 'L', 'F', 'C', 'M', 'A', 'G', 'T', 'W', 
                           'S', 'Y', 'P', 'H', 'D', 'E', 'N', 'Q', 'K', 'R'),
                group = c(4.5, 4.2, 3.8, 2.8, 2.5, 1.9, 1.8, -0.4, -0.7, -0.9, -0.8,
                          -1.3, -1.6, -3.2, -3.5, -3.5, -3.5, -3.5, -3.9, -4.5),
                stringsAsFactors=F
              ), 
              
              # Colour based on nucleotide
              nucleotide2 = data.frame(
                letter = c('A', 'C', 'G', 'T', 'U'),
                col = c('darkgreen', 'blue', 'orange', 'red', 'red'),
                stringsAsFactors = F
              ), 
              
              #alt red BA1200
              nucleotide = data.frame(
                letter = c('A', 'C', 'G', 'T', 'U'),
                col = c('#109648', '#255C99', '#F7B32B', '#D62839', '#D62839'),
                stringsAsFactors = F
              ), 
              
              base_pairing = data.frame(
                letter = c('A', 'T', 'U', 'G', 'C'),
                group = c(rep('Weak bonds', 3), rep('Strong bonds', 2)),
                col = c(rep('darkorange', 3), rep('blue', 2)),
                stringsAsFactors = F
              ),
              
              # ClustalX color scheme: 
              # http://www.jalview.org/help/html/colourSchemes/clustal.html
              clustalx = data.frame(
                letter = c('W', 'L', 'V', 'I', 'M', 'F', 'A', 'R', 'K', 'T', 'S', 'N', 'Q', 'D', 'E', 'H', 'Y', 'C', 'G', 'P'),
                col = c(rep('#197FE5', 7), rep('#E53319', 2), rep('#19CC19', 4), rep('#CC4CCC', 2), 
                        rep('#19B2B2', 2), '#E57F7F', '#E5994C', '#B0B000'),
                stringsAsFactors = F
              ),
              
              # Taylor color scheme (PMID: 9342138)
              taylor = data.frame(
                letter = c('D','S','T','G','P','C','A','V','I','L','M','F','Y','W','H','R','K','N','Q','E'),
                col = c('#FF0000','#FF3300','#FF6600','#FF9900','#FFCC00','#FFFF00','#CCFF00','#99FF00',
                        '#66FF00','#33FF00','#00FF00','#00FF66','#00FFCC','#00CCFF','#0066FF','#0000FF',
                        '#6600FF','#CC00FF','#FF00CC','#FF0066'),
                stringsAsFactors = F
              )
  )
  
  if(!'group' %in% names(cs)) cs$group = cs$letter
  
  # Set attributes
  attr(cs, 'cs_label') = col_scheme
  class(cs) = c('data.frame','ggseqlogo_cs')
  
  return(cs)
}



cs <- get_col_scheme("hydropathy")
tmp = cs[!duplicated(cs$group) & !is.na(cs$group),]
col_map = unlist( split(tmp$col, tmp$group) )


