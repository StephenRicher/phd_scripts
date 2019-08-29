library(hicrep)

data_dir = "/media/stephen/Data/hic_analysis_v2/hicexplorer/all_regions"
qc_dir = "/home/stephen/x_db/DBuck/s_richer/stephen_test/projects/hic_analysis2/qc"

capture_regions = read.csv("/home/stephen/phd/scripts/capture_regions.bed", sep = "\t", 
                           header = FALSE, col.names = c("chromosome", "start", "end", "region"))

# Initialise dataframe to store values
scc_values <- data.frame("region" = character(0), "binsize" = numeric(0),
                         "sample1" = character(0), "sample2" = character(0),
                         "h" = numeric(0), "scc" = numeric(0),
                         stringsAsFactors=FALSE)

# Define sample names with replicate number
samples = c("HB2_WT-1", "HB2_WT-2", "HB2_CL4-1", "HB2_CL4-2", "MCF7-1", "MCF7-2")

for (region in capture_regions$region) {
  
  # Set max interaction as half capture region size
  start = capture_regions[capture_regions$region == region, "start"]
  end = capture_regions[capture_regions$region == region, "end"]
  max_interaction = as.integer((end - start)/2)
  
  for (bin in rev(seq(1000, 10000, 1000))) {
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
          m2_path = paste(
            data_dir, "/", region, "/", bin, "/", sample2, "-", region, "-", bin, "_hicrep.tsv", 
            sep = "")
          # If file does not exist or exists with zero size
          if (!(file.exists(m2_path) && file.info(m2_path)$size != 0)) {
            next
          }
          m2 = read.csv(m2_path, sep = '\t', header = FALSE)
          if (sample1 == "HB2_WT-1" && sample2 == "HB2_WT-2") {
            # Calculate optimal smoothing parameter for a given region and bin
            print(
              paste('Calculating optimal smoothing parameter for region: ', region, ', bin: ', bin, '.', 
              sep = ""))
            h_hat <- htrain(m1, m2, 
                            resol = bin, 
                            max = max_interaction, 
                            range = 0:20)
          }
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
          write.csv(scc_values, paste(qc_dir, "hicrep_scc_all.csv", sep = ""))
        }
      }
    }
  }
}
write.csv(scc_values, paste(qc_dir, "hicrep_scc_all.csv", sep = ""))

