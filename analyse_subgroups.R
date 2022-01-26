






# Heatmap -----------------------------------------------------------------
subgroup_heatmap <- function(df, file) {
  mat <- subgroups_sample %>% 
    select(sample, subgroup, prop) %>% 
    spread(sample, prop, fill = 0) %>% column_to_rownames("subgroup")

  pheatmap::pheatmap(
    mat, 
    main = "",
    # filename = "results/djc_cdna/heatmap_hinge_subgroups.pdf",
    filename = file,
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
    #number_format = "%.4f",
    number_color = "black",
    fontsize_number = 8,
    fontsize = 10,
    angle_col = "90",
    cluster_rows = F,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE
  )
}




# Pie plot ----------------------------------------------------------------
subgroup_pie <- function(df, file) {
  df %>%
    mutate(label = paste0(subgroup, " (", round(prop *100, 2), "%)")) %>%
    ggplot(aes(x = 1, y = prop, fill = reorder(label, -prop))) +
    geom_bar(stat = "identity") +
    xlim(0, 2) +
    coord_polar(theta = "y") +
    theme_void() +
    theme(legend.title = element_blank()) +
    scale_fill_manual(values = colors) 
}


# ggplot(subgroups_sample, aes(x = "", y = prop, fill = subgroup)) +
#   geom_bar(stat = "identity", position = position_stack(reverse = T), width = 0.9) +
#   coord_polar(theta = "y")  +
#   facet_wrap(sample ~ ., nrow = 2) +
#   geom_text(aes(y = prop, label = label), position = position_stack(reverse = T, vjust = 0.5), size = 2.5) +
#   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
#   theme_void() +
#   theme(axis.title = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks = element_blank()
#   ) +
#   labs(x = "", y = "Proprotion", fill = "Subgroup")
# 
# p_subgroups
# 
# pdf("results/djc_cdna/hinge_subgroups.pdf")
# p_subgroups
# dev.off()