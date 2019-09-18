library(diffHic)
library(edgeR)
library(csaw)
library(rhdf5)
library(stringr)
library(statmod).
# Define directory and move to path
path="/home/stephen/x_db/DBuck/s_richer/stephen_test/projects/hic_analysis/post_alignment/diffhic/"
di_path=paste(path,"differential_interactions/", sep = "")
setwd(path)

# Whole genome sequence - FASTA headers have been modified to ensure they only contain chromosomes.
genome = paste(path, "Homo_sapiens.GRCh38.dna.primary_assembly.captured_regions.fa", sep = "")
hs.frag = cutGenome(genome, "GATC", 4)
hs.param <- pairParam(hs.frag)


capture_regions = read.csv("/home/stephen/phd/scripts/capture_regions.bed", sep = "\t", 
                           header = FALSE, col.names = c("chromosome", "start", "end", "region"))

hicPlotTADs = "/home/stephen/miniconda3/bin/hicPlotTADs"
track_template_path = "/home/stephen/phd/scripts/pyGenomeTracks_configs/"
HB2_WT_vs_HB2_CL4_template_tads.ini = paste(track_template_path, "diffhic_HB2_WTvsHB2_CL4_tads.ini", sep = "")
HB2_WT_vs_MCF7_template_tads.ini = paste(track_template_path, "diffhic_HB2_WTvsMCF7_tads.ini", sep = "")
HB2_WT_vs_HB2_CL4_template.ini = paste(track_template_path, "diffhic_HB2_WTvsHB2_CL4.ini", sep = "")
HB2_WT_vs_MCF7_template.ini = paste(track_template_path, "diffhic_HB2_WTvsMCF7.ini", sep = "")


# Samples to process.
samples = c("HB2_WT_1", "HB2_WT_2", "HB2_CL4_1", "HB2_CL4_2", "MCF7_1", "MCF7_2")

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

bin.size = 5000
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
mval <- adj.counts[,1]-adj.counts[,5]
smoothScatter(ab, mval, xlab="A", ylab="M", main="KO (1) vs. Flox (2)")
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")     

## Modelling Biological Variablity ##

# Coef 2 measures HB2_WT vs HB2_CL4 and coef 3 measures HB2_WT vs MCF7
design <- model.matrix(
  ~factor(c("HB2_WT","HB2_WT","HB2_CL4", "HB2_CL4", 
            "MCF7", "MCF7"),
          levels = c("HB2_WT", "HB2_CL4", "MCF7")))
colnames(design) <- c("intercept", "HB2_CL4", "MCF7")

# Convert to DGEList for edgeR
y <- asDGEList(data)
y <- estimateDisp(y, design)
plotBCV(y)
plotMDS(y)
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)

for (coef in c(2, 3)) {
  if(coef == 2) {
    comp = "HB2_CL4"
  } else {
    comp = "MCF7"
  }
  result <- glmQLFTest(fit, coef = coef)
  rowData(data) <- cbind(rowData(data), result$table)
  
  clustered.sig <- diClusters(data, result$table, target=0.01, cluster.args=list(tol=1))
  useful.cols <- as.vector(outer(c("seqnames", "start", "end"), 1:2, paste0))
  tabcom <- combineTests(clustered.sig$indices[[1]], result$table)
  tabbest <- getBestTest(clustered.sig$indices[[1]], result$table)
  tabstats <- data.frame(tabcom[,1:4], logFC=tabbest$logFC, FDR=clustered.sig$FDR)
  results.d <- as.data.frame(clustered.sig$interactions)[,useful.cols]
  results.d <- cbind(results.d, tabstats)
  o.d <- order(results.d$PValue)
  write.table(results.d[o.d,], file = paste(di_path, "HB2_WT-vs-", comp, ".tsv", sep = ""), sep="\t",quote=FALSE, row.names=FALSE)
  
  for (tads in c('tads', 'logFC')) {
    if (tads == 'tads') {
      if (comp == "HB2_CL4") {
        pyGenome_template.ini = HB2_WT_vs_HB2_CL4_template_tads.ini 
      } else {
        pyGenome_template.ini = HB2_WT_vs_HB2_CL4_template_tads.ini
      }
    } else {
      if (comp == "HB2_CL4") {
        pyGenome_template.ini = HB2_WT_vs_HB2_CL4_template.ini 
      } else {
        pyGenome_template.ini = HB2_WT_vs_HB2_CL4_template.ini
      }
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
      
        pyGenome.ini = paste(di_path, "HB2_WT_vs_", comp, "-", region, "-", bin.size, "-", tads,".ini", sep = "")
        
        system(paste("sed 's/capture_region/", region, "/g; s/binsize/", bin.size, "/g' ", pyGenome_template.ini, " > ", pyGenome.ini, sep = ""))
      
        # Remove section for config file if no significant UP/DOWN interactions detected.
        if (nrow(arc[arc$logFC < -2,]) != 0) {
          write.table(arc[arc$logFC < -2,], file = paste(di_path, region, "_HB2_WT-vs-", comp, "_di_down.arc", sep = ""), sep = "\t",quote = FALSE, 
                      row.names = FALSE, col.names = FALSE)
        } else {
          system(paste("sed -i '/Start DI_Down/,/End DI_Down/d' ", pyGenome.ini, sep = ""))
        }
        if (nrow(arc[arc$logFC > 2,]) != 0) {
          write.table(arc[arc$logFC > 2,], file = paste(di_path, region, "_HB2_WT-vs-", comp, "_di_up.arc", sep = ""), sep = "\t",quote = FALSE, 
                      row.names = FALSE, col.names = FALSE)
        } else {
          system(paste("sed -i '/Start DI_Up/,/End DI_Up/d' ", pyGenome.ini, sep = ""))
          # Replace line below # DI_Down share axis with height config to prevent overlap with HiC map.
          system(paste("sed -i '/# DI_Down share axis/{n;s/.*/height = 10/}' ", pyGenome.ini, sep = ""))
        }
        # Run hicPlotTADs only if atleast 1 DI interaction detected
        if (nrow(arc[abs(arc$logFC) > 2,]) != 0) {
          system(paste(hicPlotTADs, " --tracks ", pyGenome.ini, " --region ", chromosome, ":", start, "-", end, " -o ", di_path, region, "_HB2_WT_vs_", comp, "-", tads, ".png", sep = ""))
        }
      }
    }
  }
}


















