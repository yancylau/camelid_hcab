library(ape)
library(ggtree)
library(ggmsa)
library(Biostrings)
library(ConsensusClusterPlus)

load("data/djc_cdna/vdjc.Rda")


ch2 <- vdjc %>% 
  subset(stop_codon == FALSE) %>%
  subset(!is.na(ch2_seq)) %>%
  group_by(ch2_seq) %>% summarise(count = n()) %>%
  arrange(-count) %>%
  rownames_to_column("seq_id") %>%
  mutate(seq_id = paste0("ch2_", seq_id)) %>%
  mutate(ch2_seq = str_remove_all(ch2_seq, "-")) %>%
  subset(nchar(ch2_seq) == 70) %>% 
  mutate(prop = count / sum(count)) 


ggplot(ch2, aes(x = seq_id, y = prop)) +
  geom_bar(stat = "identity") 





x <- AAStringSet(ch2$ch2_seq) 
names(x) <- ch2$seq_id

d <- as.dist(stringDist(x)/width(x)[1])

d2 <- as.matrix(d)
cons_clu_data = as.matrix(d)

cluster_result = ConsensusClusterPlus(
  cons_clu_data,
  maxK = 9,
  reps = 1000,
  pItem = 0.8,
  pFeature = 1,
  title = "results/djc_cdna/cnsensusclustering_ch2",
  clusterAlg = "km",
  distance = "euclidean",
  # clusterAlg = "hc",
  #distance = "pearson",
  seed = 123,
  plot = "png"
)


# Select cluster
n_cluster <- 5

consensus_class <- as.matrix(cluster_result[[n_cluster]][["consensusClass"]])
sample_order <- cluster_result[[n_cluster]][["consensusTree"]]$order
clustered_samples <- consensus_class[sample_order, ]
table(clustered_samples)
 


ordered_d <- d2[sample_order,]

anno_col <- data.frame(clustered_samples) %>%
  rownames_to_column("sample") %>%
  dplyr::rename(cluster = clustered_samples) %>%
  mutate(cluster = paste0("C", cluster)) %>%
  column_to_rownames("sample")

anno_colors = list(Cluster = c(C1 = "red", C2 = "blue", C3 = "green", C4 = "yellow"))
pheatmap(
  filename = "results/djc_cdna/heatmap_of_ch2_clustered.pdf",
  1-ordered_d,
  cluster_row = T,
  cluster_cols = T,
  show_colnames = T,
  show_rownames = F,
  annotation_colors = anno_colors,
  annotation_col = anno_col,
  # # annotation_row = annotation_row,
  # fontsize = 6.5,
  # fontsize_row = 6,
  fontsize_col = 2,
  # gaps_col = 50
  angle_col = "90",
)



reorderd_sample_names <- rownames(ordered_d)
ch2 %>%
  mutate(seq_id = factor(seq_id, levels = reorderd_sample_names )) %>%
  ggplot(aes(x = seq_id, y = prop)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, colour = "black", hjust=0.95,vjust=0.3),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(color = "black"))



fviz_nbclust(d2, kmeans, method = "wss")

km <- kmeans(d2, centers = 2, nstart = 25)

fviz_cluster(
  km, 
  data = d2, 
  ellipse.type = "norm", 
  pointsize = .5,
  geom = c("point")
)

library("FactoMineR")
res.pca <- PCA(d2, graph = FALSE)



dat = d2
dat$pc1 <- res.pca$ind$coord[,1] # indexing the first column
dat$pc2 <- res.pca$ind$coord[,2]# indexing the first column


ggplot(dat, aes(x = pc1, y = pc2)) +
  geom_point()


fviz_pca_var(res.pca, col.var = "black")

iris.pca <- PCA(iris[,-5], graph = FALSE)

fviz_pca_ind(iris.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = iris$Species, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)






data("diamonds")

# put data into a dataframe (rather than a tibble)
dat <- diamonds %>% data.frame

dat <- dat[sample(rownames(dat), 2000),]


dat <-  dat %>% filter (x > 0, y > 0, z > 0)

# log price

# center and scale the data
for (i in 1:length(colnames(dat))){
  
  if (is.numeric(dat[, i])==TRUE)
    
    dat[, i] <- as.numeric(scale(dat[, i]))
  
  else
    
    dat[, i] <- dat[, i]
  
}


pca1 <- PCA(dat[ ,c("carat", "depth", "table", "price", "x", "y", "z", "clarity", "cut", "color")], 
            
            quali.sup = c(8:10), graph = FALSE)


plot.PCA(pca1)


dat$pc1 <- pca1$ind$coord[, 1] # indexing the first column

dat$pc2 <- pca1$ind$coord[, 2]  #

pca.vars <- pca1$var$coord %>% data.frame

pca.vars$vars <- rownames(pca.vars)

pca.vars.m <- melt(pca.vars, id.vars = "vars")

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

circ <- circleFun(c(0,0),2,npoints = 500)


p <- ggplot(data = dat, aes(x = pc1, y = pc2)) +
  
  geom_point()


levels(dat$cut)
## [1] "Fair"      "Good"      "Very Good" "Premium"   "Ideal"
p <- ggplot(data = dat, aes(x = pc1, y = pc2, color = cut, shape = cut)) +
  
  geom_hline(yintercept = 0, lty = 2) +
  
  geom_vline(xintercept = 0, lty = 2) +
  
  guides(color = guide_legend(title = "Cut"), shape = guide_legend(title = "Cut")) +
  
  scale_shape_manual(values = c(15, 16, 16, 17, 18)) +
  
  geom_point(alpha = 0.8, size = 2) 


p

p <- p + stat_ellipse(geom="polygon", aes(fill = cut), 
                      
                      alpha = 0.2, 
                      
                      show.legend = FALSE, 
                      
                      level = 0.95) +
  
  xlab("PC 1 (68.25%)") + 
  
  ylab("PC 2 (18.37%)") +
  
  theme_minimal() +
  
  theme(panel.grid = element_blank(), 
        
        panel.border = element_rect(fill= "transparent"))

p
