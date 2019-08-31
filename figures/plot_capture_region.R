#!/usr/bin/Rscript

library("ggplot2")
library("RColorBrewer")
library("ggpubr")

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("Two argument must be supplied: input data and output directory", call. = FALSE)
}
dir = args[2]
if (!endsWith(dir, '/')) {
  dir = paste(dir, '/',  sep = '')
}
dir.create(paste(dir, 'capture_region_density', sep = ''))
width = 10
height = (9/16) * width

capture_regions = read.csv(args[1], sep = "\t", header = FALSE, 
                           col.names = c("chromosome", "start", "end", "region"))

data = data.frame()
bin_range = seq(1, 50, 1)
for (region in capture_regions$region) {
  start = capture_regions[capture_regions$region == region, 'start']
  end = capture_regions[capture_regions$region == region, 'end']
  length = end - start
  nbins = (length / (bin_range*1000))
  tmp_data = data.frame(bin_range = bin_range, nbins = nbins, region = region)
  data = rbind.data.frame(data, tmp_data)
}

colourCount = length(unique(capture_regions$region))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

g = ggplot(data, aes(x = bin_range, y = nbins, colour = region)) + 
  geom_line() +
  scale_y_log10() +
  scale_colour_manual(name = NULL,
                      values = colorRampPalette(brewer.pal(8, "Dark2"))(colourCount)) +
  labs(title = NULL, subtitle = NULL, tag = NULL,
       x = "Bin size (kb)", y = "Number of bins", caption = NULL) +
  theme_pubr(legend = "right", base_size = 12) +
  guides(colour = guide_legend(ncol = 1))
ggsave(filename = paste(dir,'capture_region_density/nbins_by_binsize.png', sep = ''), 
       plot = g, 
       dpi = 300, width = width, height = height)
