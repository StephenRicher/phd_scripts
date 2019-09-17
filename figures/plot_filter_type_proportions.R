#!/usr/bin/Rscript

library("ggplot2")
library("RColorBrewer")
library("ggpubr")
library(dplyr)
library(tidyr)
library(stringr)

filter_types = read.csv('/media/stephen/Data/hic_analysis/data/filter_statistics.tsv',
                        sep = '\t', header = FALSE)

filter_types = cbind.data.frame(filter_types, 
                                str_split_fixed(filter_types$V1, "-", 2))

names(filter_types) = c('full_sample_name', 'Filter', 'Count', 'Sample', 'Replicate')

# Reorder levels by number because level name is variable i.e < 1000bp
filter_types$Filter = factor(filter_types$Filter,
                             levels(filter_types$Filter)[c(9,3,4,8,1,2,5,6,7)])

# Remove total read count
filter_types = filter_types[filter_types$Filter != "Total" & filter_types$Filter != "Filtered" & filter_types$Count != 0 ,]

# Set order for 'Sample' column
#filter_types$Sample = factor(filter_types$Sample, levels = c('HB2 WT','HB2 TSS KO','MCF7'))

width = 10
height = (9/16) * width
dir = '~/phd/phd_report/figures/'
# Expand colour palette to number of groups
ggplot(filter_types, 
       aes(x = Replicate, y = Count,
           fill = Filter)) + 
  geom_bar(position = 'fill', stat = 'identity') +
  facet_wrap(~Sample) +
  scale_fill_brewer(name = NULL, palette = "Dark2", direction = 1) + 
  labs(title = NULL, subtitle = NULL, tag = NULL,
       x = 'Replicate', y = 'Percentage of read pairs',
       caption = NULL) +
  theme_pubr(legend = 'right', base_size = 14) +
  guides(colour = guide_legend(nrow = 1)) +
  scale_y_continuous(labels = scales::percent_format()) + 
  ggsave(filename = paste(dir,"filter_type_proportion.png", sep = ""), 
         dpi = 300, width = width, height = height)

