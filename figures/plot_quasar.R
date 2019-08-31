#!/usr/bin/Rscript
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("Two argument must be supplied: input data and output directory", call. = FALSE)
}
dir = args[2]
if (!endsWith(dir, '/')) {
  dir = paste(dir, '/',  sep = '')
}
dir.create(paste(dir, 'quasar_qc', sep = ''))
width = 10
height = (9/16) * width

data = read.table(args[1], col.names = c('sample', 'binsize', 'scale', 'coverage', 'score', 'score2', 'region'))

data = data[data$score != 'NaN',]
data$binsize = as.numeric(as.character(data$binsize))
data$score = as.numeric(as.character(data$score))
data = data.frame(lapply(data, function(x) {gsub("HB2_CL4", "HB2_TSS-KO", x)}))

for (region in unique(data$region)) {
  quasar = ggplot(data[data$region == region,], 
                  aes(x = as.numeric(as.character(binsize)), 
                      y = as.numeric(as.character(score)), 
                      colour = sample)) +
    geom_line() + geom_point() +
    scale_colour_brewer(name = NULL, palette = "Dark2") + 
    labs(title = NULL, subtitle = NULL, tag = NULL,
         x = "Bin size (kb)", y = "Quasar QC score", caption = NULL) +
    theme_pubr(legend = "bottom") +
    guides(colour = guide_legend(nrow = 1))
  ggsave(filename = paste(dir,'quasar_qc/quasar-', region, '.png', sep = ''), 
         plot = quasar, 
         dpi = 300, width = width, height = height)
}