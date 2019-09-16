#!/usr/bin/Rscript

library("ggplot2")
library("RColorBrewer")
library("ggpubr")
library(dplyr)
library(tidyr)

truncate_stats = read.csv('/media/stephen/Data/hic_analysis_v2/qc/hictools-truncate_summary.txt', header = FALSE)
names(truncate_stats) = c('sample', 'statistic', 'count')

# Set order for 'Sample' column
truncate_stats$sample = factor(truncate_stats$sample, levels = rev(c('HB2_TSS_KO-1-R1', 'HB2_TSS_KO-1-R4',
                                                                 'HB2_TSS_KO-2-R1', 'HB2_TSS_KO-2-R4',
                                                                 'HB2_WT-1-R1', 'HB2_WT-1-R4',
                                                                 'HB2_WT-2-R1', 'HB2_WT-2-R4',
                                                                 'MCF7-1-R1', 'MCF7-1-R4',
                                                                 'MCF7-2-R1', 'MCF7-2-R4')))

# Remove total count
truncate_stats = truncate_stats[truncate_stats$statistic != 'Total',]
truncate_stats = truncate_stats[truncate_stats$statistic != 'Mean truncated length',]


width = 10
height = (9/16) * width
dir = '~/phd/phd_report/figures/'
# Expand colour palette to number of groups
ggplot(truncate_stats, 
       aes(x = sample, y = count,
           fill = statistic)) + 
  geom_bar(position = 'fill', stat = 'identity') +
  scale_fill_brewer(name = NULL, palette = "Dark2") +
  labs(title = NULL, subtitle = NULL, tag = NULL,
       x = NULL, y = NULL,
       caption = NULL) +
  theme_pubr(legend = "bottom") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(labels = scales::percent_format()) +
  coord_flip() +
  ggsave(filename = paste(dir,"truncation_stats.png", sep = ""), 
         dpi = 300, width = width, height = height)
