#!/usr/bin/env Rscript


source("R/variability_index.R")


load("data/v_cdna/v_cdna.Rda")
load("data/v_gdna/v_gdna.Rda")


# Germline V genes

# Rearranged V genes


# Calculate variability index table
indexes_v_gdna <- make_indexes_table(v_gdna)
indexes_v_cdna <- make_indexes_table(v_cdna)

# Plot variability
variability_plot_gdna <- plot_variability(indexes_v_gdna)
variability_plot_cdna <- plot_variability(indexes_v_cdna)

pdf("results/camelids/variability_plot_classical.pdf", width = 6, height = 4)
variability_plot_gdnas
variability_plot_cdna
dev.off()









# # Save
# pdf("results/camelids/variability_plot.pdf", width = 6, height = 4)
# p_variability
# dev.off()


# pdf("results/trend3.pdf", width=12, height=6)
# plot_variability(v_indexes)
# plot_variability_change(v_indexes)
# plot_variability_diff(v_indexes)
# dev.off()
# 
