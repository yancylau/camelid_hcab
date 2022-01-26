#!/usr/bin/env Rscript
#
# 1. Filter sequences by sub-region length 
# 2. Select fields
# 3. IMGT numbering for nt and aa sequences
#
# Raw (n=60,659) -> functional (n=60,653) -> region_length (n=56,011)


dir.create("data/v_cdna", showWarnings = F)

# Alignment results to the germline Vs of Bactrian camel
igblast_bactrian <- read_tsv("imported/v_cdna/igblast/bactrian/atleast3.airr") 
selected_fields = c("sequence_id", "sample", "species", "dupcount", "sequence", "stop_codon", 
                    "v_call", "v_identity", "sequence_alignment", "sequence_alignment_aa", 
                    "fwr1", "cdr1", "fwr2", "cdr2", "fwr3",
                    "fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa", "fwr3_aa")
igblast <- igblast_bactrian %>% 
  mutate(sequence_id = str_remove_all(sequence_id, "SAMPLE=|DUPCOUNT=")) %>% 
  separate(sequence_id, c("sequence_id", "sample", "dupcount"), sep = ";") %>%
  mutate(sample = recode(sample, 
                         B1 = "cDNA-dome-1",
                         B4 = "cDNA-dome-2",
                         B5 = "cDNA-dome-3",
                         B6 = "cDNA-dome-4",
                         B7 = "cDNA-dome-5",
                         Y3 = "cDNA-wild-1",
                         Y4 = "cDNA-wild-2",
                         Y5 = "cDNA-wild-3")) %>%
  mutate(species = ifelse(sample %in% c("cDNA-wild-1", "cDNA-wild-2", "cDNA-wild-3"), "wild", "dome")) %>%
  select(all_of(selected_fields)) %>%
  separate(v_call, into = "v_call", sep = ",", extra = "drop")
  

  
# Plot region length
p_region_length <- igblast %>% 
  mutate(
    fr1_length = nchar(fwr1_aa),
    fr2_length = nchar(fwr2_aa),
    fr3_length = nchar(fwr3_aa),
    cdr1_length = nchar(cdr1_aa),
    cdr2_length = nchar(cdr2_aa),
  ) %>%
  select(sample, fr1_length, fr2_length, fr3_length, cdr1_length, cdr2_length) %>%
  gather(key = "region", value = "length", -sample) %>%
  group_by(sample, region, length) %>% summarise(count = n()) %>%
  group_by(sample, region) %>% mutate(prop = count / sum(count)) %>%
  ggplot(aes(x = length, y = prop)) +
  geom_bar(stat = "identity") +
  facet_grid(sample ~ region, scales = "free", space = "free") + 
  theme_void()
p_region_length

pdf("results/v_cdna/region_length_pre-filter.pdf")
p_region_length
dev.off()


#### Filter by region length and IMGT numbering 
source("R/imgt_numbering.R")

igblast <- igblast %>% 
  subset(stop_codon == F) %>%
  subset(nchar(fwr1_aa) %in% c(25, 26)) %>%
  subset(nchar(fwr2_aa) %in% c(17)) %>%
  subset(nchar(fwr3_aa) %in% c(38, 39)) %>%
  subset(nchar(cdr1_aa) < 13) %>%
  subset(nchar(cdr2_aa) < 11) %>%
  mutate(
    fwr1_imgt_aa = fr_numbering_aa(fwr1_aa, "FR1"),
    fwr2_imgt_aa = fr_numbering_aa(fwr2_aa, "FR2"),
    fwr3_imgt_aa = fr_numbering_aa(fwr3_aa, "FR3"),
    cdr1_imgt_aa = cdr_numbering_aa(cdr1_aa, "CDR1"),
    cdr2_imgt_aa = cdr_numbering_aa(cdr2_aa, "CDR2")
  ) %>%
  mutate(
    fwr1_imgt = fr_numbering_nt(fwr1, "FR1"),
    fwr2_imgt = fr_numbering_nt(fwr2, "FR2"),
    fwr3_imgt = fr_numbering_nt(fwr3, "FR3"),
    cdr1_imgt = cdr_numbering_nt(cdr1, "CDR1"),
    cdr2_imgt = cdr_numbering_nt(cdr2, "CDR2")
  ) %>%
  unite("imgt_aa", c(fwr1_imgt_aa, cdr1_imgt_aa, fwr2_imgt_aa, cdr2_imgt_aa, fwr3_imgt_aa), remove = FALSE, sep = "") %>%
  unite("imgt_nt", c(fwr1_imgt, cdr1_imgt, fwr2_imgt, cdr2_imgt, fwr3_imgt), remove = FALSE, sep = "") 



save(igblast, file = "data/v_cdna/igblast.Rda")




# http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
#
# IMGT unique numbering (aa/nt)
# FWR1 :     1 -  26 (26)      1-78
# CDR1 :    27 -  38 (12)     81-114
# FWR2 :    39 -  55 (17)    117-165
# CDR2 :    56 -  65 (10)    168-195
# FWR3 :    66 - 104 (39)    198-312
# CDR3 :   105 - 117 (13)    313-351
# FWR4 :   118 - 129（12）   352-287


