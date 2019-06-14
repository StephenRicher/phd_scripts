library(diffHic)
library(edgeR)
library(csaw)
library(rhdf5)
library(stringr)

capture_regions = read.csv("/home/stephen/h/phd/scripts2/hic_scripts/capture_regions.bed", sep = "\t", 
                           header = FALSE, col.names = c("chromosome", "start", "end", "region"))

# Define capture region to analysis
region="GNG12_AS1_DIRAS3"

chromosome = capture_regions[capture_regions$region == region, "chromosome"]
start = capture_regions[capture_regions$region == region, "start"]
end = capture_regions[capture_regions$region == region, "end"]
location = paste(chromosome, start+1, end, sep = "-")

setwd(paste("/home/stephen/x_db/DBuck/s_richer/stephen_test/projects/hic_analysis/post_alignment/diffhic/", region, sep = ""))

# Generate restriction digest
genome = paste("Homo_sapiens.GRCh38.dna.primary_assembly.", region, "-", location, ".fa", sep = "")
hs.frag = cutGenome(genome, "GATC", 4)
hs.param <- pairParam(hs.frag)

samples = c("HB2_WT_1", "HB2_WT_2", "HB2_CL4_1", "HB2_CL4_2", "MCF7_1", "MCF7_2")
for (sample in samples) {
  bam = paste(sample, ".", region, ".bam", sep = "")
  matrix = paste(sample, ".h5", sep = "")
  matrix_trim = paste(sample, "-trim.h5", sep = "")
  sample_diagnostics = preparePairs(bam, hs.param, file = matrix, 
                                    dedup = TRUE, minq = 20)
  prunePairs(matrix, hs.param, file.out = matrix_trim)
  print(paste(sample,sample_diagnostics))
}

bin.size = 5000

data <- squareCounts(paste(samples, "-trim.h5", sep = ""), hs.param, filter = 0, width = bin.size)

## FILTERING UNINTERESTING INTERACTIONS ##

# Remove low abundance counts
ave.ab <- aveLogCPM(asDGEList(data))
hist(ave.ab, xlab="Average abundance", col="grey80", main="")
count.keep = ave.ab >= aveLogCPM(5, lib.size=mean(data$totals))
summary(count.keep)

# Remove diagonal
ndiag.keep <- filterDiag(data)
summary(ndiag.keep)

margin.data <- marginCounts(input, hs.param, width = bin.size)
#margin.data

adjc <- cpm(asDGEList(margin.data), log = TRUE, prior.count = 1)
smoothScatter(0.5*(adjc[,3]+adjc[,4]), adjc[,3]-adjc[,4],xlab="A", ylab="M", main="HB2_WT(1) vs. MCF7(1)")

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
mval <- adj.counts[,3] - adj.counts[,4]
smoothScatter(ab, mval, xlab="A", ylab="M", main="KO (1) vs. Flox (2)")
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")

data <- normOffsets(data, method = "loess", se.out=TRUE)
nb.off <- assay(data, "offset")

ab <- aveLogCPM(asDGEList(data))
o <- order(ab)
adj.counts <- log2(assay(data) + 0.5) - nb.off/log(2)
mval <- adj.counts[,1]-adj.counts[,3]
smoothScatter(ab, mval, xlab="A", ylab="M", main="KO (1) vs. Flox (2)")
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")     

## Modelling Biological Variablity ##

# Coef 2 measures HB2_WT vs HB2_CL4 and coef 3 measures HB2_WT vs MCF7
design <- model.matrix(~factor(c("HB2_WT","HB2_WT","HB2_CL4", "HB2_CL4", "MCF7", "MCF7"),
                               levels = c("HB2_WT", "HB2_CL4", "MCF7")))
colnames(design) <- c("intercept", "HB2_CL4", "MCF7")



# Convert to DGEList for edgeR
y <- asDGEList(data)
y <- estimateDisp(y, design)
plotBCV(y)
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
  write.table(results.d[o.d,], file = paste(region, "_HB2_WT-vs-", comp, ".tsv", sep = ""), sep="\t",quote=FALSE, row.names=FALSE)
  
  results.d[,c("start2", "end2","start1", "end1")] = results.d[,c("start2", "end2","start1", "end1")] + start
  arc=results.d[,c("seqnames2","start2", "end2", "seqnames1","start1", "end1", "logFC")]
  
  write.table(arc[arc$logFC < 0 & abs(arc$logFC) > 2,], file = paste(region, "_HB2_WT-vs-", comp, "_di_down.arc", sep = ""), sep = "\t",quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
  write.table(arc[arc$logFC > 0 & abs(arc$logFC) > 2,], file = paste(region, "_HB2_WT-vs-", comp, "_di_up.arc", sep = ""), sep = "\t",quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
}


#hicPlotTADs --tracks hic_track_GNG12-AS1-WT_vs_CL4.ini --region 1:65834316-70234314 -o GNG12-AS1_HB2_WT_vs_HB2_CL4.png
#hicPlotTADs --tracks hic_track_GNG12-AS1-WT_vs_MCF7.ini --region 1:65834316-70234314 -o GNG12-AS1_HB2_WT_vs_MCF7.png



#hicPlotTADs --tracks hic_track_TET1-WT_vs_CL4.ini --region 10:66560360-70694482 -o TET1_HB2_WT_vs_HB2_CL4.png
#hicPlotTADs --tracks hic_track_TET1_vs_MCF7.ini --region 10:66560360-70694482 -o TET1_HB2_WT_vs_MCF7.png

























x=abs(results.d$logFC)
results.d$score = as.integer((1000-0)*((x-min(x))/(max(x)-min(x))) + 0)
a = data.frame(paste("chr", results.d[o.d, c("seqnames1")], sep = ""),
               do.call(pmin, results.d[o.d, c("start1", "start2")]),
               do.call(pmax, results.d[o.d, c("end1", "end2")]),
               ".", 
               results.d[o.d, c("score")],
               results.d[o.d, c("logFC")],
               ".", "0",
               paste("chr", results.d[o.d, c("seqnames1")], sep = ""),
               do.call(pmin, results.d[o.d, c("start1", "start2")]),
               do.call(pmin, results.d[o.d, c("end1", "end2")]),
               ".", ".",
               paste("chr", results.d[o.d, c("seqnames2")], sep = ""),
               do.call(pmax, results.d[o.d, c("start1", "start2")]),
               do.call(pmax, results.d[o.d, c("end1", "end2")]),
               ".", ".")
colnames(a) = c("track type=interact name=\"interact Example Two\" description=\"Chromatin interactions\" useScore=on maxHeightPixels=50:100:200 visibility=full",
                rep("", ncol(a)-1))
write.table(a, file = "HB2_WT_vs_HB2_CL4.tsv", sep = "\t",quote = FALSE, 
            row.names = FALSE, col.names = TRUE)
