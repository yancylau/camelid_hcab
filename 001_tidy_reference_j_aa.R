
dir.create("data/refs", showWarnings = F)


ighj <- read_csv("doc/refs/IMGT_GENE_DB/IGHJ/fr4.csv")

# Drop sequences whose FR4 is not complete
j <- ighj %>% subset(complete==T) %>%
  mutate(gene = str_extract(id, "IGHJ[0-9]+")) %>%
  select(species, id, gene, aa, fr4, complete) %>%
  mutate(molecule = "gDNA")



save(j, file = "data/refs/j.Rda")




# 
# 
# 
# 
# ighj_web <- read_csv("../refs/IMGT_GENE_DB/IGHJ/ighj_web.csv") %>% rename_all(tolower) 
# 
# j_web <- ighj_web %>% 
#   subset(functionality == "F") %>% 
#   subset(complete == T)
# 
# 
# j_web %>% group_by(species) %>% summarise(count = n())
# 
# 
# j_web %>% 
#   mutate(a = str_sub(fr4, 6, 6)) %>%
#   mutate(gene = str_extract(gene, "IGHJ[0-9]+")) %>%
#   group_by(species, gene, a) %>% summarise(count = n()) -> d
