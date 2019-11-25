#!/usr/bin/Rscript

library(hicrep)
library(pheatmap)
library(reshape2)
library(RColorBrewer)
library(igraph)
library(ggplot2)
library(viridisLite)

args = commandArgs(trailingOnly=TRUE)

if (args[1] == 'allele') {
  samples = c("HB2_WT_G1-1", "HB2_WT_G1-2", "HB2_WT_G2-1", "HB2_WT_G2-2")
  qc_dir = '/home/stephen/x_db/DBuck/s_richer/hic_01/allele_specific/qc'
  data_dir = '/home/stephen/x_db/DBuck/s_richer/hic_01/allele_specific/hicexplorer/all_regions'
  binrange = seq(5000, 20000, 1000) 
} else {
  samples = c("HB2_WT-1", "HB2_WT-2", "HB2_TSS-KO 1", "HB2_TSS-KO-2", "MCF7-1", "MCF7-2")
  qc_dir = '/home/stephen/x_db/DBuck/s_richer/hic_01/qc/'
  data_dir = '/home/stephen/x_db/DBuck/s_richer/hic_01/data/hicexplorer/all_regions/'
  binrange = seq(1000, 10000, 1000)
}

# All hicrep output will be stored in a folder in qc directory
dir.create(paste(qc_dir, 'hicrep', sep = '/'))
capture_regions = read.csv("/home/stephen/phd/scripts/capture_regions.bed", sep = "\t", 
                           header = FALSE, col.names = c("chromosome", "start", "end", "region"))

# Initialise dataframe to store values
scc_values <- data.frame("region" = character(0), "binsize" = numeric(0),
                         "sample1" = character(0), "sample2" = character(0),
                         "h" = numeric(0), "scc" = numeric(0),
                         stringsAsFactors=FALSE)

for (region in capture_regions$region) {
  
  # Set max interaction as half capture region size, or 1,000,000bp
  start = capture_regions[capture_regions$region == region, "start"]
  end = capture_regions[capture_regions$region == region, "end"]
  max_interaction = min(1000000, as.integer((end - start)/2))
  
  for (bin in binrange) {
    h_hat = NULL
    sample_combinations = c()
    for (sample1 in samples) {
      m1_path = paste(
        data_dir, "/", region, "/", bin, "/", sample1, "-", region, "-", bin, "_hicrep.tsv", 
        sep = "")
      # If file does not exist or exists with zero size
      if (!(file.exists(m1_path) && file.info(m1_path)$size != 0)) {
        next
      }
      m1 = read.csv(m1_path, sep = '\t', header = FALSE)
      for (sample2 in samples) {
        if (sample1 != sample2) {
          #next
          m2_path = paste(
            data_dir, "/", region, "/", bin, "/", sample2, "-", region, "-", bin, "_hicrep.tsv", 
            sep = "")
          # If file does not exist or exists with zero size
          if (!(file.exists(m2_path) && file.info(m2_path)$size != 0)) {
            next
          }
          m2 = read.csv(m2_path, sep = '\t', header = FALSE)
          # Calculate optimal smoothing parameter for a given region and bin
          print(
            paste('Calculating optimal smoothing parameter for region: ', region, ', bin: ', bin, '.', sep = ""))
            tryCatch({h_hat <- htrain(m1, m2, 
                                      resol = bin, 
                                      max = max_interaction, 
                                      range = 0:20)},
                     error = function(e) {cat("ERROR :",conditionMessage(e), "\n"); h_hat = NULL})
            
          
          # h_hat will be NULL if HB2_WT-1 or HB2_WT-2 dont exists for region and bin
          if (is.null(h_hat)) {
            next
          }
          
          # Don't rerun equivalent combinations
          combo = paste(sort(c(sample1, sample2)), collapse = ":")
          if (combo %in% sample_combinations){
            next
          } else {
            sample_combinations = c(sample_combinations, combo)
          }
          
          print(paste(
            'Calculating SCC for region: ', region, ', bin: ', bin, ', against ', sample1, ' and ', sample2, '.', 
            sep = ""))
          
          # Downsample matrices to the same sequencing depth
          min_sample = min(sum(m1[,-c(1:3)]), sum(m2[,-c(1:3)]))
          m1_ds <- depth.adj(m1, min_sample, bin, out = 0) 
          m2_ds <- depth.adj(m2, min_sample, bin, out = 0) 
        
          # Format HiC matrix pairs and smooth
          pre_hic <- prep(m1_ds, m2_ds, 
                          resol = bin, 
                          h = h_hat, 
                          max = max_interaction)
          
          # Calculate SCC
          SCC.out = get.scc(pre_hic, 
                            resol = bin, 
                            max = max_interaction)
          
          print(paste(
            'SCC for region: ', region, ', bin: ', bin, ', against ', sample1, ' and ', sample2, ' is ', SCC.out$scc,
            sep = ""))
          
          scc_values[nrow(scc_values) + 1,] = list(region, bin, sample1, sample2, h_hat, SCC.out$scc)
          write.csv(scc_values, paste(qc_dir, 'hicrep/hicrep_scc_all.csv', sep = '/'))
        } else {
          scc_values[nrow(scc_values) + 1,] = list(region, bin, sample1, sample2, NA, NA)
        }
      }
    }
    subset = scc_values[scc_values$region == region & scc_values$binsize == bin, 
                        c('sample1', 'sample2', 'scc')]
    if (nrow(subset) == 0) {
      next
    }
    g <- graph.data.frame(subset, directed=FALSE)
    heatmap = pheatmap(get.adjacency(g, attr="scc", sparse=FALSE),
                       color = viridis(100), display_numbers = TRUE,
                       number_color = 'red', 
                       fontsize = 8, fontsize_number = 10,
                       angle_col = 0, treeheight_col = 0)
    ggsave(paste(qc_dir, '/hicrep/hicrep-', region, '-', bin, '.png', sep = ''), 
           plot = heatmap, width = 8, height = 5)
  }
}
write.csv(scc_values, paste(qc_dir, 'hicrep/hicrep_scc_all.csv', sep = '/'))

for (region in capture_regions$region) {
  for (bin in binrange) {
    subset = scc_values[scc_values$region == region & scc_values$binsize == bin, 
                        c('sample1', 'sample2', 'scc')]
    if (nrow(subset) == 0) {
      next
    }
    g <- graph.data.frame(subset, directed=FALSE)
    heatmap = pheatmap(get.adjacency(g, attr="scc", sparse=FALSE),
                       color = viridis(100), display_numbers = TRUE,
                       number_color = 'red', 
                       fontsize = 8, fontsize_number = 10,
                       angle_col = 0, treeheight_col = 0)
    ggsave(paste(qc_dir, '/hicrep/hicrep-', region, '-', bin, '.png', sep = ''), 
           plot = heatmap, width = 8, height = 5)
  }
}
