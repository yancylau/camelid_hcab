

# Representative genes (OTUs) distribution across samples


load("data/v_gdna/v_gdna.Rda")


clusters <- read_tsv("imported/v_gdna/igblast/atleast3_0.99/clustering/clusters.tsv") %>% 
  mutate(otu = str_remove_all(otu, "otu"))

types <- v_gdna %>% select(sequence_id, otu, type)

otus <- clusters %>% 
  left_join(types, by = "otu") %>%
  group_by(sample, sequence_id, otu, type) %>% 
  summarise(size = n(), count = sum(count)) %>%
  ungroup() %>%
  # mutate(label = paste0(size, "(", count, ")")) %>%
  mutate(sample = recode(sample, 
                         GZ2 = "gDNA-dome-1",
                         SZ4 = "gDNA-dome-2",
                         JTR1 = "gDNA-dome-3",
                         JTE2 = "gDNA-dome-4",
                         JTE3 = "gDNA-dome-5",
                         JTE4 = "gDNA-dome-6",
                         JTE5 = "gDNA-dome-7",
                         YTE2 = "gDNA-wild-1")) 

save(otus, file = "data/v_gdna/otus.Rda")
write_csv(otus, "results/v_gdna/otus.csv")




#### Upset R
d <- otus %>%
  mutate(size = 1) %>%
  select(-c(otu, count, type)) %>%
  spread(sample, size, fill = 0) %>%
  #mutate(otu = paste0("ighv", otu)) %>%
  column_to_rownames("sequence_id")
p_upset <- upset(d, 
                 nsets = 8, 
                 # sets = c("JTE3", "JTE2", "JTR1", "SZ4", "GZ2", "JTE4", "JTE5", "YTE2"),
                 # sets = c("YTE2", "JTE5", "JTE4", "GZ2", "SZ4", "JTR1", "JTE2", "JTE3"),
                 keep.order = TRUE,
                 order.by = "degree",  
                 sets.x.label = "Number of detected genes",
                 mainbar.y.label = "Common genes")

p_upset
pdf("results/v_gdna/germline_v_distribution_upsetr.pdf",  width = 5, height = 4)
p_upset 
dev.off()



# 
# #### PCA analysis
# library(FactoMineR)
# library(factoextra)
# 
# d <- otus %>%
#   select(-c(otu, type, count)) %>%
#   spread(sample, size, fill = 0) %>%
#   column_to_rownames("sequence_id")
# res.pca <- PCA(d, graph = FALSE)
# 
# fviz_pca_ind(res.pca, col.ind = "cos2", pointsize = "cos2",
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE # Avoid text overlapping (slow if many points)
# )
# 
# 
# col_ind <- c("domestic", "domestic", "domestic", "domestic", "domestic", "domestic","domestic", "wild")
# fviz_pca_ind(res.pca,
#              # geom.ind = "point", # show points only (nbut not "text")
#              # col.ind = col_ind, # color by groups
#              # palette = c("#00AFBB", "#E7B800"),
#              addEllipses = TRUE, # Concentration ellipses
#              legend.title = "Groups"
# )
# 
# 
# #### Sequence length
# #### BLAST result (identity)
# #### PCA analysis
# #### 饱和曲线、稀疏曲???






