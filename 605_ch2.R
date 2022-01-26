library(ape)
library(ggtree)
library(ggmsa)
library(Biostrings)

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

## Different aa > 3, prop > 0.1
ch2_subgroups <- vdjc %>% 
  subset(stop_codon == FALSE) %>%
  subset(!is.na(ch2_seq)) %>%
  group_by(ch2_seq) %>% summarise(count = n()) %>%
  arrange(-count) %>% 
  head(2) %>% select(ch2_seq) %>%
  mutate(ch2_subgroup = c("ch2_1", "ch2_2"))

save(ch2_subgroups, file = "data/djc_cdna/ch2_subgroups.Rda")


d <- vdjc %>% 
  left_join(ch2_subgroups, by = "ch2_seq") %>% subset(!is.na(ch2_subgroup)) %>%
  group_by(sample, ch2_subgroup) %>% summarise(count = n()) %>%
  group_by(sample) %>% mutate(prop = count / sum(count)) %>%
  arrange(sample, -prop) %>%
  mutate(label = paste0(round(prop, 3) * 100, "%"))



p_subgroups <- ggplot(d, aes(x = "", y = prop, fill = ch2_subgroup)) +
  geom_bar(stat = "identity", position = position_stack(reverse = T), width = 0.9) +
  coord_polar(theta = "y")  +
  facet_wrap(sample ~ ., nrow = 2) +
  geom_text(aes(y = prop, label = label), position = position_stack(reverse = T, vjust = 0.5), size = 2.5) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
  theme_void() +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()
  ) +
  labs(x = "", y = "Proprotion", fill = "Subgroup")

p_subgroups

pdf("results/djc_cdna/ch2_subgroups.pdf")
p_subgroups
dev.off()





### heatmap of hinge subgrouo 
d2 <- d %>% select(-count, -label) %>%
  spread(sample, prop, fill = 0) %>%
  column_to_rownames("ch2_subgroup")


pheatmap::pheatmap(
  d2, 
  main = "",
  filename = "results/djc_cdna/heatmap_ch2_subgroups.pdf",
  width = 8, height = 8,
  # cellheight = 4,
  # cellweight = 0.05,
  na_col = "white",
  # breaks = seq(0.8, 1, 0.05),
  # scale = "none",
  # # color = colorRampPalette(rev(brewer.pal(n = 7, name)))
  # # legend_breaks = seq(0.8, 1, 0.02),
  # annotation_row = anno_row,
  # annotation_col = anno_col,
  # # annotation_colors = anno_colors,
  display_numbers = TRUE,
  # number_format = "%.4f",
  number_color = "black",
  fontsize_number = 8,
  fontsize = 10,
  angle_col = "90",
  cluster_rows = FALSE,
  cluster_cols = F,
  show_rownames = TRUE,
  show_colnames = TRUE
)


