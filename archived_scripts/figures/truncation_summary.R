library("ggplot2")
library("RColorBrewer")
library("plotly")
library("stringr")
library("ggpubr")

setwd("/home/stephen/x_db/DBuck/s_richer/stephen_test/projects/hic_analysis/qc_logs/")


truncation_summary = read.table("all_samples.trunc_summary.txt",
                                 col.names = c('sample', 'length', 'count'))

ggplot(truncation_summary,
       aes(x = length, y = count, colour = sample)) +
  # Normalize each facet panel seperately by total count of that facet group
  geom_density(stat = 'identity') + 
  scale_colour_brewer(name = NULL, palette = "Dark2") + 
  labs(title = NULL, subtitle = NULL, tag = NULL,
       x = "Insert size (bp)", y = "Frequency", caption = NULL) +
  theme_pubr(legend = "bottom")
