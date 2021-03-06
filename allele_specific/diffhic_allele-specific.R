library(diffHic)
library(edgeR)
library(csaw)
library(rhdf5)
library(stringr)
library(statmod)

# Define directory and sub folders and move to path
path="/home/stephen/x_db/DBuck/s_richer/hic_01/allele_specific/diffhic/"
di_path = paste(path,"differential_interactions/", sep = "")
tracks_path = paste(path,"tracks/", sep = "")
figs_path = paste(path,"figures/", sep = "")
dir.create(di_path)
dir.create(tracks_path)
dir.create(figs_path)
setwd(path)

# Whole genome sequence - FASTA headers have been modified to ensure they only contain chromosomes.
genome = paste(path, "Homo_sapiens.GRCh38.dna.primary_assembly.captured_regions.fa", sep = "")
hs.frag = cutGenome(genome, "GATC", 4)
hs.param <- pairParam(hs.frag)

capture_regions = read.csv("/home/stephen/phd/scripts/capture_regions.bed", sep = "\t", 
                           header = FALSE, col.names = c("chromosome", "start", "end", "region"))

hicPlotTADs = "/home/stephen/anaconda3/envs/hic_analysis/bin/hicPlotTADs"
track_template_path = "/home/stephen/phd/scripts/pyGenomeTracks_configs/"
HB2_WT_G1_vs_G2_template.ini = paste(track_template_path, "diffhic_G1vsG2.ini", sep = "")
HB2_WT_G1_vs_G2_template_tads.ini = paste(track_template_path, "diffhic_G1vsG2_tads.ini", sep = "")

# Samples to process.
samples = c("HB2_WT_G1-1",  "HB2_WT_G2-1", "HB2_WT_G1-2", "HB2_WT_G2-2") 

# Generate HiC matrix of trans interactions within capture regions.
for (sample in samples) {
  print(sample)
  bam = paste(sample, ".captured.bam", sep = "")
  matrix = paste(sample, ".h5", sep = "")
  matrix_trim = paste(sample, ".trim.h5", sep = "")
  if (!file.exists(matrix_trim)) {
    sample_diagnostics = preparePairs(bam, hs.param, file = matrix, 
                                      dedup = TRUE, minq = 20)
    prunePairs(matrix, hs.param, file.out = matrix_trim)
    print(paste(sample,sample_diagnostics))
  }
}

# Thresholds for filtering diffential interactions
logfc_threshold = 1
fdr_threshold = 0.05

bin.size = 10000
input = paste(samples, ".trim.h5", sep = "")
data <- squareCounts(input, hs.param, width = bin.size)
colnames(data) = samples
margin.data <- marginCounts(input, hs.param, width=bin.size)
adjc <- cpm(asDGEList(margin.data), log=TRUE, prior.count=5)
colnames(adjc) = samples
smoothScatter(0.5*(adjc[,1]+adjc[,2]), adjc[,1]-adjc[,2],xlab="A", ylab="M", main="Flox (1) vs. Ko (1)")

# Remove low abundance counts
ave.ab <- aveLogCPM(asDGEList(data))
hist(ave.ab, xlab="Average abundance", col="grey80", main="")
count.keep = ave.ab >= aveLogCPM(5, lib.size=mean(data$totals))
summary(count.keep)

# Remove diagonal
ndiag.keep <- filterDiag(data)
summary(ndiag.keep)

trended <- filterTrended(data)
smoothScatter(trended$log.distance, trended$abundances,xlab="Log-Distance", ylab="Normalized abundance")
o <- order(trended$log.distance)
lines(trended$log.distance[o], trended$threshold[o], col="red", lwd=2)
trend.keep <- trended$abundances > trended$threshold
summary(trend.keep)

original <- data
data <- data[trend.keep & ndiag.keep & count.keep,]

## Normalisation ##

#cnv.offs <- normalizeCNV(data, margin.data)
#assay(data, "offset") <- cnv.offs# updating the offset matrix

ab <- aveLogCPM(asDGEList(data))
o <- order(ab)
adj.counts <- cpm(asDGEList(data), log = TRUE)
mval <- adj.counts[,1] - adj.counts[,3]
smoothScatter(ab, mval, xlab="A", ylab="M", main="KO (1) vs. Flox (2)")
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")

data <- normOffsets(data, method = "loess", se.out=TRUE)
nb.off <- assay(data, "offset")

ab <- aveLogCPM(asDGEList(data))
o <- order(ab)
adj.counts <- log2(assay(data) + 0.5) - nb.off/log(2)
mval <- adj.counts[,1]-adj.counts[,4]
smoothScatter(ab, mval, xlab="A", ylab="M", main="KO (1) vs. Flox (2)")
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")     

## Modelling Biological Variablity ##

# Coef 2 measures HB2_WT vs HB2_CL4 and coef 3 measures HB2_WT vs MCF7

## Modelling Biological Variablity ##
sample <- factor(c("HB2_WT-1", "HB2_WT-1", "HB2_WT-2", "HB2_WT-2"))
genome <- factor(c("G1", "G2", "G1", "G2"), levels = c("G1", "G2"))
design <- model.matrix(~sample+genome)

# Convert to DGEList for edgeR
y <- asDGEList(data)
y <- estimateDisp(y, design)
plotBCV(y)
plotMDS(y)
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)

result <- glmQLFTest(fit)
rowData(data) <- cbind(rowData(data), result$table)
  
clustered.sig <- diClusters(data, result$table, target = fdr_threshold, cluster.args=list(tol=1))
useful.cols <- as.vector(outer(c("seqnames", "start", "end"), 1:2, paste0))
tabcom <- combineTests(clustered.sig$indices[[1]], result$table)
tabbest <- getBestTest(clustered.sig$indices[[1]], result$table)
tabstats <- data.frame(tabcom[,1:4], logFC=tabbest$logFC, FDR=clustered.sig$FDR)
results.d <- as.data.frame(clustered.sig$interactions)[,useful.cols]
results.d <- cbind(results.d, tabstats)
o.d <- order(results.d$PValue)
write.table(results.d[o.d,], file = paste(di_path, "HB2_WT_G1-vs-G2.tsv", sep = ""), sep="\t",quote=FALSE, row.names=FALSE)
  
for (tads in c('logFC', 'tads')) {
  if (tads == 'tads') {
    pyGenome_template.ini = HB2_WT_G1_vs_G2_template_tads.ini
  } else {
    pyGenome_template.ini = HB2_WT_G1_vs_G2_template.ini
  }
    
  for (region in capture_regions$region) {
    chromosome = capture_regions[capture_regions$region == region, "chromosome"]
    start = capture_regions[capture_regions$region == region, "start"]
    end = capture_regions[capture_regions$region == region, "end"]
    location = paste(chromosome, start+1, end, sep = "-")
      
    region_results.d = results.d[results.d$seqnames1 == region,]
      
    if (nrow(region_results.d) != 0) {
      region_results.d[,c("start2", "end2","start1", "end1")] = region_results.d[,c("start2", "end2","start1", "end1")] + start
      region_results.d[,c("seqnames1", "seqnames2")] = chromosome
      arc=region_results.d[,c("seqnames2","start2", "end2", "seqnames1","start1", "end1", "logFC")]
        
      pyGenome.ini = paste(tracks_path, "HB2_WT_G1-vs-G2", "-", region, "-", bin.size, "-", tads,".ini", sep = "")
        
      system(paste("sed 's/capture_region/", region, "/g; ",
                        "s/binsize/", bin.size, "/g; ", 
                        "s/logfc_threshold/", logfc_threshold, "/g; ",
                        "s/fdr_threshold/", fdr_threshold, "/g' ", pyGenome_template.ini, " > ", pyGenome.ini, sep = ""))
        
      # Remove section for config file if no significant UP/DOWN interactions detected.
      if (nrow(arc[arc$logFC < -logfc_threshold,]) != 0) {
        write.table(arc[arc$logFC < -logfc_threshold,], file = paste(di_path, region, "_HB2_WT_G1-vs-G2_di_down.arc", sep = ""), sep = "\t",quote = FALSE, 
                    row.names = FALSE, col.names = FALSE)
      } else {
        system(paste("sed -i '/Start DI_Down/,/End DI_Down/d' ", pyGenome.ini, sep = ""))
      }
      if (nrow(arc[arc$logFC > logfc_threshold,]) != 0) {
        write.table(arc[arc$logFC > logfc_threshold,], file = paste(di_path, region, "_HB2_WT_G1-vs-G2_di_up.arc", sep = ""), sep = "\t",quote = FALSE, 
                    row.names = FALSE, col.names = FALSE)
      } else {
        system(paste("sed -i '/Start DI_Up/,/End DI_Up/d' ", pyGenome.ini, sep = ""))
        # Replace line below # DI_Down share axis with height config to prevent overlap with HiC map.
        system(paste("sed -i '/# DI_Down share axis/{n;s/.*/height = 10/}' ", pyGenome.ini, sep = ""))
        }
      # Run hicPlotTADs only if atleast 1 DI interaction detected
      if (nrow(arc[abs(arc$logFC) > logfc_threshold,]) != 0) {
        system(paste(hicPlotTADs, " --tracks ", pyGenome.ini, " --region ", chromosome, ":", start, "-", end, " -o ", figs_path, region, "HB2_WT_G1-vs-G2", "-", tads, ".png", sep = ""))
      }
    }
  }
}
















