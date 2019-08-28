#!/usr/bin/Rscript

library("ggplot2")
library("RColorBrewer")
library("ggpubr")
library(dplyr)
library(tidyr)

filter_stats = read.csv('/media/stephen/Data/hic_analysis_v2/filter_stats.csv')

# Remove total read count
filter_stats = filter_stats[filter_stats$Filter != "Total reads",]

filter_stats$Replicate = factor(filter_stats$Replicate)

# Set order for 'Sample' column
filter_stats$Sample = factor(filter_stats$Sample, levels = c('HB2 WT','HB2 TSS KO','MCF7'))

filter_stats$Filter = factor(filter_stats$Filter, 
                             levels = c('Reads <20bp','Both unmapped','Only R1 mapped',
                                        'Only R2 mapped', 'Duplicates', 'MAPQ <15',
                                        'HiC filtered', 'Not both captured', 'Retained'))

width = 10
height = (9/16) * width
dir = '~/phd/phd_report/figures/'
# Expand colour palette to number of groups
colours <- colorRampPalette(brewer.pal(8, "Dark2"))(9)
ggplot(filter_stats, 
       aes(x = Replicate, y = Count,
           fill = Filter)) + 
  geom_bar(position = 'fill', stat = 'identity') +
  facet_wrap(~Sample) +
  #scale_fill_manual(values = brewer.pal(9,"Greys"),name = NULL) +
  scale_fill_manual(values = colours, name = NULL) +
  labs(title = NULL, subtitle = NULL, tag = NULL,
       x = "Replicate", y = "Number of read pairs",
       caption = NULL) +
  theme_pubr(legend = "right") +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(labels = scales::percent_format()) +
  ggsave(filename = paste(dir,"filter_stats.png", sep = ""), 
         dpi = 300, width = width, height = height)
