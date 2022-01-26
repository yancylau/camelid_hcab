
fasta <- read_tsv("../imported/refs/imgt_ighv.tab", col_names = F) %>% 
  select(1:2) %>% set_names("id", "imgt_nt")

ndm <- read_tsv("../imported/refs/human.ndm.imgt", col_names = F) %>%
  set_names(c("id", "fwr1_start", "fwr1_end", 
              "cdr1_start", "cdr1_end", 
              "fwr2_start", "fwr2_end", 
              "cdr2_start", "cdr2_end", 
              "fwr3_start", "fwr3_end", 
              "chain", "stop")) %>%
  mutate(type = "nt")

pdm <- read_tsv("../imported/refs/human.pdm.imgt", col_names = F) %>%
  set_names(c("id", "fwr1_start", "fwr1_end", 
              "cdr1_start", "cdr1_end", 
              "fwr2_start", "fwr2_end", 
              "cdr2_start", "cdr2_end", 
              "fwr3_start", "fwr3_end", 
              "chain", "stop")) %>%
  mutate(type = "aa")


dm <- bind_rows(ndm, pdm) %>% subset(str_detect(id, "IGHV|VH")) %>%
  group_by(id) %>% add_count(name = "count") 

dm %>% subset(count != 2) -> d0
  

d <- full_join(dm, fasta) %>%
  mutate(length = nchar(imgt_nt))

d %>% subset(count != 2 & length == 0) -> d







