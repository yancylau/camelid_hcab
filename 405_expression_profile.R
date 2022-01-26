#!/usr/bin/env Rscript


load("data/v_cdna/v_cdna.Rda")
load("data/v_gdna/v_gdna.Rda")


v_cdna = v_cdna %>% subset(type != "unknown")

get_c_positions <- function(v) {
  position_list <- v %>% pull(imgt_aa) %>% str_locate_all("C")
  
  lapply(seq_along(position_list), function(i) {
    id = v$id[i]
    pos = position_list[[i]][,1]
    data.frame(id = id, position = pos)
  }) %>% bind_rows() 
}

c_groups <- get_c_positions(v_gdna) %>% 
  group_by(id) %>% add_count(name = "count") %>%
  mutate(flag = 1) %>%
  spread(position, flag, fill = 0) %>% 
  unite(group, matches("[0-9]+"), sep = "") 



## Alignment identity
v_identity <- v_cdna %>% 
  select(type, v_identity, v_call) %>%
  separate(v_call, into = "v_call", sep = ",", extra = "drop")

p_identity <- ggplot(v_identity, aes(x =  type, y = v_identity, color = type)) +
  geom_violin(aes(fill = type)) +
  geom_boxplot(width=0.1) + 
  # geom_jitter(width = 0.1, size = 0.1) +
  theme_classic() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 0, colour = "black", hjust=0.95, vjust=0.3),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "", y = "Identity", color = "Type")
p_identity 


## Heatmap
count_matrix <- v_cdna %>% 
  # subset(type == "VHH") %>%
  separate(v_call, into = "v_call", sep = ",", extra = "drop") %>%
  group_by(sample, v_call) %>% summarise(count = n()) %>%
  spread(sample, count, fill = 0) %>%
  column_to_rownames("v_call") 

anno_row <- v_cdna %>%  select(v_call) %>%
  separate(v_call, into = "v_call", sep = ",", extra = "drop") %>%  distinct() %>%
  separate(v_call, into = c(NA, "type"), sep = "_", extra = "drop", remove = F) %>%
  left_join(c_groups, by = c("v_call" = "id")) %>%
  select(-count) %>%
  column_to_rownames("v_call")
  
anno_col <- v_cdna %>% select(sample, species) %>% distinct() %>% column_to_rownames("sample")
anno_colors <- list(type = c(vhh = "red", vh = "blue"))

pheatmap::pheatmap(
  log2(count_matrix + 1), 
  filename = "results/v_cdna/heatmap_of_v_calls.pdf",
  main = "",
  # cellheight = 4,
  # cellweight = 0.002,
  na_col = "white",
  border_color = NA,
  # breaks = seq(0.8, 1, 0.05),
  # scale = "none",
  # scale = "row",
  #scale = "column",
  #legend_breaks = seq(0.8, 1, 0.02),
  annotation_row = anno_row,
  annotation_col = anno_col,
  annotation_colors = anno_colors,
  angle_col = "90",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  display_numbers = FALSE
)



