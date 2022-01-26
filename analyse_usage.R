# Heatmap, Pie plot, Barplot to show gene usage


source("R/colors.R")





# Heatmap -----------------------------------------------------------------
usage_heatmap <- function(df, file) {
  anno_col <- data.frame(
    sample = c("cDNA-dome-1", "cDNA-dome-2", "cDNA-dome-3", "cDNA-dome-4",
               "cDNA-dome-5", "cDNA-wild-1", "cDNA-wild-3"), 
    species = c(rep("domestic", 5), rep("wild", 2))
  ) %>%
    column_to_rownames("sample")
  
  mat <- df %>% 
    ungroup() %>% select(sample, subgene, prop) %>% 
    spread(sample, prop, fill = 0) %>% 
    column_to_rownames("subgene")
    
  pheatmap::pheatmap(
    mat, 
    main = "",
    filename = file,
    width = 8, height = 8,
    # cellheight = 4,
    # cellweight = 0.05,
    na_col = "white",
    # breaks = seq(0.8, 1, 0.05),
    # scale = "row",
    # # color = colorRampPalette(rev(brewer.pal(n = 7, name)))
    # # legend_breaks = seq(0.8, 1, 0.02),
    # annotation_row = anno_row,
    annotation_col = anno_col,
    # # annotation_colors = anno_colors,
    display_numbers = TRUE,
    number_format = "%.4f",
    number_color = "black",
    fontsize_number = 8,
    fontsize = 10,
    angle_col = "90",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE
  )
}




# Pie plot ----------------------------------------------------------------
usage_pie <- function(df, file) {
  df %>%
    mutate(label = paste0(subgene, " (", round(prop *100, 2), "%)")) %>%
    ggplot(aes(x = 1, y = prop, fill = reorder(label, -prop))) +
    geom_bar(stat = "identity") +
    xlim(0, 2) +
    coord_polar(theta = "y") +
    theme_void() +
    theme(legend.title = element_blank()) +
    scale_fill_manual(values = colors) 
}




# Bar plot ----------------------------------------------------------------
usage_bar <- function(df, file) {
  ggplot(df, aes(x = subgene, y = mean, fill = subgene)) +
    geom_bar(stat="identity") +
    geom_errorbar(aes(x=subgene, ymin=mean-sd, ymax=mean+sd), width = .5) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2))
}






# Seqlogo -----------------------------------------------------------------

# Color scheme
hydropathy_col =  ggseqlogo::make_col_scheme(
  chars = c(
    c("A", "C", "I", "L", "M", "F", "W", "V"), 
    c("R", "N", "D", "Q", "E", "K"),
    c("G", "H", "P", "S", "T", "Y")
  ),
  groups = c(rep("Hydrophobic", 8), rep("Hydrophilic", 6), rep("Neutral", 6)),
  cols = c(rep("#109648", 8), rep("#F7B32B", 6), rep("#D62839", 6))
)

seqlogo_plot <- function(seqs) {
  ggplot() + 
    ggseqlogo::geom_logo(seqs, method = "prob", col_scheme = hydropathy_col) + 
    ggseqlogo::theme_logo() + 
    theme(legend.position = "right") +
    # theme(axis.text.y = element_text(size = 6)) +
    NULL 
}



