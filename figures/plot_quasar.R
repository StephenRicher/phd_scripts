#!/usr/bin/Rscript
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop('Three argument must be supplied: quasar output, captured regions and output directory', 
       call. = FALSE)
}

# Create quasar qc folder in output dir
outdir = args[3]
if (!endsWith(outdir, '/')) {
  outdir = paste(outdir, '/',  sep = '')
}
dir.create(paste(outdir, 'quasar_qc', sep = ''))

# Define dimensions for each plot
width = 10
height = (9/16) * width

# Read quasar QC data
quasar_data = read.table(args[1], 
                         col.names = c('sample', 'binsize', 'scale', 'coverage', 'score', 'score2', 'region'))
quasar_data = quasar_data[quasar_data$score != 'NaN',]
quasar_data$binsize = as.numeric(as.character(quasar_data$binsize))
quasar_data$score = as.numeric(as.character(quasar_data$score))
quasar_data = data.frame(lapply(quasar_data, function(x) {gsub("HB2_CL4", "HB2_TSS-KO", x)}))

# Read capture regions and calculate nbins per length
capture_regions = read.csv(args[2], sep = "\t", header = FALSE, 
                           col.names = c("chromosome", "start", "end", "region"))
capture_stats = data.frame()
max_bin = 50
bin_range = seq(1, max_bin, 1)
for (region in capture_regions$region) {
  start = capture_regions[capture_regions$region == region, 'start']
  end = capture_regions[capture_regions$region == region, 'end']
  length = end - start
  nbins = (length / (bin_range*1000))
  tmp_data = data.frame(bin_range = bin_range, nbins = nbins, region = region)
  capture_stats = rbind.data.frame(capture_stats, tmp_data)
}


# Plot graphs for each region
for (region in unique(quasar_data$region)) {
  quasar = ggplot(quasar_data[quasar_data$region == region,], 
                  aes(x = as.numeric(as.character(binsize)), 
                      y = as.numeric(as.character(score)), 
                      colour = sample)) +
    geom_line() + geom_point() +
    scale_colour_brewer(name = NULL, palette = "Dark2") + 
    labs(title = NULL, subtitle = NULL, tag = 'a)',
         x = "Bin size (kb)", y = "Quasar QC score", caption = NULL) +
    theme_pubr(legend = "bottom", base_size = 14) +
    guides(colour = guide_legend(nrow = 1))
  
  capture = ggplot(data = capture_stats[capture_stats$region == region, ], 
                   aes(x = bin_range, y = nbins, colour = region)) + 
    geom_line() +
    scale_y_log10() +
    labs(title = NULL, subtitle = NULL, tag = 'b)',
         x = "Bin size (kb)", y = "Number of bins", caption = NULL) +
    theme_pubr(base_size = 14) +
    theme(legend.position = "none")
  
  ggsave(filename = paste(outdir,'quasar_qc/quasar-', region, '.png', sep = ''), 
         grid.arrange(rbind(ggplotGrob(quasar), 
                            ggplotGrob(capture))), 
         device = "png", dpi = 300, width = width, height = height * 2)
}




