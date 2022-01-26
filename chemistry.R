#!/usr/bin/env Rscript
#


## IMGT 'Physicochemical' classes of the 20 common amino acids 
# http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html
# J Mol Recognit. Jan-Feb 2004;17(1):17-32. doi: 10.1002/jmr.647.


# 'Hydropathy'
hydrophobic = c("A", "C", "I", "L", "M", "F", "W", "V")
neutral = c("G", "H", "P", "S", "T", "Y")
hydrophilic = c("R", "N", "D", "Q", "E", "K")
hydropathy <- data.frame(
  composite = c(hydrophobic,  neutral, hydrophilic),
  hydropathy = c(rep("Hydrophobic", 8), rep("Neutral", 6), rep("Hydrophilic", 6))
)

# 'Volume'
very_sampll = c("A", "G", "S")
small= c("N", "D", "C", "P", "T")
medium = c("Q", "E", "H", "V")
large = c("R", "I", "L", "K", "M")
very_large =  c("F", "W", "Y")
volume <- data.frame(
  composite = c(very_sampll, small, medium, large, very_large),
  # volume = c(rep(1, 3), rep(2, 5), rep(3, 4), rep(4, 5), rep(5, 3),
  volume = c(rep("very_samll", 3), rep("small", 5), rep("medium", 4), rep("large", "5"), rep("very_large", 3))
)

# 'Chemical'
aliphatic = c("A", "G", "I", "L", "P", "V")
aromatic = c("F", "W", "Y")
sulfur = c("C", "M")
hydroxyl = c("S", "T")
basic = c("R", "H", "K")
acidic = c("D", "E")
amide = c("N", "Q")
chemical <- data.frame(
  composite = c(aliphatic, aromatic, sulfur, hydroxyl, basic, acidic, amide),
  chemical = c(rep("aliphatic", 6), rep("aromatic", 3), rep("sulfur", 2), rep("hydroxyl", 2), 
               rep("basic", 3), rep("acidic", 2), rep("amide", 2))
)

# 'Charge'
positive = c("R", "H", "K")
negative = c("D", "E")
uncharged = c("A", "N", "C", "Q", "G", "I", "L", "M", "F", "P", "S", "T", "W", "Y", "V")
charge <- data.frame(
  composite = c(positive, negative, uncharged),
  charge = c(rep("positive", 3), rep("negative", 2), rep("uncharged", 15))
)

# 'Hydrogen donor or acceptor atoms'
donor = c("R", "K", "W")
acceptor = c("D", "E")
donor_and_acceptor = c("N", "Q", "H", "S", "T", "Y")
none = c("A", "C", "G", "I", "L", "M", "F", "P", "V")
hydrogen <- data.frame(
  composite = c(donor, acceptor, donor_and_acceptor, none),
  hydrogen = c(rep("donor", 3), rep("acceptor", 2), rep("donor_and_acceptor", 6), rep("none", 9))
)

# 'Polarity'
polar = c("R", "N", "D", "Q", "E", "H", "K", "S", "T", "Y")
nonpolar = c("A", "C", "G", "I", "L", "M", "F", "P", "W", "V")
polarity <- data.frame(
  composite = c(polar, nonpolar),
  polarity = c(rep("polar", 10), rep("nonpolar", 10))
)

## Merged physicochemical
chemistry <- hydropathy %>% left_join(volume) %>% left_join(chemical) %>% 
  left_join(charge) %>% left_join(hydrogen) %>% left_join(polarity) 
