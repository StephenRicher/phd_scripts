library("ggplot2")
library("RColorBrewer")
library("plotly")
library("stringr")
library("ggpubr")
setwd("/media/stephen/Data/hic_analysis/diffhic/")
setwd("/home/stephen/x_db/DBuck/s_richer/stephen_test/projects/hic_analysis/post_alignment")
data = read.csv("all_samples.hic_info_sampled.txt", header = TRUE, sep = '\t')

data$replicate = gsub("^.*_", "", data$sample)
data$sample_group = substr(data$sample, 1 , nchar(as.character(data$sample)) - 2)

trans_stats = table(as.vector(data[data$interaction_type == "trans", "orientation"]), 
                    as.vector(data[data$interaction_type == "trans", "sample"]))
sweep(trans_stats,2,colSums(trans_stats),`/`)

ggplot(data[data$interaction_type == "cis", ],
       aes(x = insert_size, colour = orientation)) +
  # Normalize each facet panel seperately by total count of that facet group
  geom_freqpoly(aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]), bins = 100) + 
  geom_vline(xintercept = 1000, linetype = "dotted", color = "#1B9E77") +
  facet_grid(sample_group ~ replicate, scales = "free_y") +
  scale_x_continuous(trans = 'log10') +
  scale_colour_brewer(name = NULL, palette = "Dark2") + 
  labs(title = NULL, subtitle = NULL, tag = NULL,
       x = "Insert size (bp)", y = "Frequency", caption = NULL) +
  theme_pubr(legend = "bottom")
#g <- ggplot_build(i)
#unique(g$data[[1]]["colour"])

gg <- ggplotly()
gg


ggplot(data[, c("ditag_length", "sample")], 
           aes(x = ditag_length, colour = sample)) + 
  geom_density() + 
  geom_vline(xintercept = 1000, linetype = "dotted", color = "#1B9E77") +
  scale_x_continuous(trans = 'log10') +
  scale_colour_brewer(name = NULL, palette = "Dark2") + 
  labs(title = NULL, subtitle = NULL, tag = NULL,
       x = "Ditag length (bp)", y = "Frequency",
       caption = "Ditag size distribution of 10M random read pairs.") +
  theme_pubr(legend = "bottom")
dd_int <- ggplotly(dynamicTicks=T)
dd_int
