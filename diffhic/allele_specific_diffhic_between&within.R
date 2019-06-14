library(diffHic)
library(edgeR)
library(csaw)
library(stringr)
region="GNG12_AS1_DIRAS3"
location="1-65834317-70234314"
# 0-based region start
region_start = 65834316

setwd(paste("/media/stephen/Data/hic_analysis/diffhic_allele_spec/diffhic/", region, sep = ""))

# Generate restriction digest
genome = paste("Homo_sapiens.GRCh38.dna.primary_assembly.", region, "-", location, ".fa", sep = "")
hs.frag = cutGenome(genome, "GATC", 4)
hs.param <- pairParam(hs.frag)

samples = c("HB2_WT_G1-1",  "HB2_WT_G2-1", "HB2_WT_G1-2", "HB2_WT_G2-2", 
            "HB2_CL4_G1-1", "HB2_CL4_G2-1", "HB2_CL4_G1-2", "HB2_CL4_G2-2") 
for (sample in samples) {
  bam = paste(sample, ".", region, ".bam", sep = "")
  matrix = paste(sample, ".h5", sep = "")
  matrix_trim = paste(sample, "-trim.h5", sep = "")
  sample_diagnostics = preparePairs(bam, hs.param, file = matrix, 
                                    dedup = TRUE, minq = 20)
  print(paste(sample, sample_diagnostics$pairs["mapped"]), sep = '\t')
  #prunePairs(matrix, hs.param, file.out = matrix_trim)
}
#25000 works!
bin.size = 5000
data <- squareCounts(paste(samples, "-trim.h5", sep = ""), 
                     hs.param, filter = 0, width = bin.size)

## FILTERING UNINTERESTING INTERACTIONS ##

# Remove low abundance counts
ave.ab <- aveLogCPM(asDGEList(data))
hist(ave.ab, xlab="Average abundance", col="grey80", main="")
count.keep = ave.ab >= aveLogCPM(3, lib.size=mean(data$totals))
summary(count.keep)

# Remove diagonal
ndiag.keep <- filterDiag(data)
summary(ndiag.keep)

margin.data <- marginCounts(input, hs.param, width = bin.size)
margin.data


adjc <- cpm(asDGEList(margin.data), log = TRUE, prior.count = 1)
smoothScatter(0.5*(adjc[,3]+adjc[,4]), adjc[,3]-adjc[,4],xlab="A", ylab="M", main="HB2_WT(1) vs. MCF7(1)")

## BEST NOT TO INCLUDE TREND.KEEP FOR ALLELE?

trended <- filterTrended(data)
smoothScatter(trended$log.distance, trended$abundances,xlab="Log-Distance", ylab="Normalized abundance")
o <- order(trended$log.distance)
lines(trended$log.distance[o], trended$threshold[o], col="red", lwd=2)
trend.keep <- trended$abundances > trended$threshold
summary(trend.keep)

original <- data
data <- data[trend.keep & ndiag.keep & count.keep,]

## Normalisation ##



ab <- aveLogCPM(asDGEList(data))
o <- order(ab)
adj.counts <- cpm(asDGEList(data), log = TRUE)
mval <- adj.counts[,1] - adj.counts[,4]
smoothScatter(ab, mval, xlab="A", ylab="M", main="KO (1) vs. Flox (2)")
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")

data <- normOffsets(data, method = "loess", se.out=TRUE)
nb.off <- assay(data, "offset")

cnv.offs <- normalizeCNV(data, margin.data)
assay(data, "offset") <- cnv.offs# updating the offset matrix

ab <- aveLogCPM(asDGEList(data))
o <- order(ab)
adj.counts <- log2(assay(data) + 0.5) - nb.off/log(2)
mval <- adj.counts[,1]-adj.counts[,3]
smoothScatter(ab, mval, xlab="A", ylab="M", main="KO (1) vs. Flox (2)")
fit <- loessFit(x=ab, y=mval)
lines(ab[o], fit$fitted[o], col="red")     

## Modelling Biological Variablity ##
cell_lines = str_remove(samples, "_G[12]-[12]")
sample <- gl(length(unique(cell_lines)),
             2, length = length(cell_lines))
cell <- factor(cell_lines)
genome <- factor(rep(c("G1","G2"), length(cell_lines)/2))
data.frame(Cell,Sample,Treat)
design <- model.matrix(~cell+cell:sample+cell:genome)

# Convert to DGEList for edgeR
y <- asDGEList(data)
y <- estimateDisp(y, design)
plotBCV(y)
plotMDS(y)
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)

# Select coef 5:6 or just 5
result <- glmQLFTest(fit, coef="cellHB2_WT:genomeG2")
topTags(result)

# Genes that respond differentently between HB2_WT and HB2_CL4
result2 = glmQLFTest(fit, contrast=c(0,0,0,0,-1,1))
topTags(result2)

rowData(data) <- cbind(rowData(data), result$table)
clustered.sig <- diClusters(data, result$table, target=0.2, cluster.args=list(tol=1))
useful.cols <- as.vector(outer(c("seqnames", "start", "end"), 1:2, paste0))
tabcom <- combineTests(clustered.sig$indices[[1]], result$table)
tabbest <- getBestTest(clustered.sig$indices[[1]], result$table)
tabstats <- data.frame(tabcom[,1:4], logFC=tabbest$logFC, FDR=clustered.sig$FDR)
results.d <- as.data.frame(clustered.sig$interactions)[,useful.cols]
results.d <- cbind(results.d, tabstats)
o.d <- order(results.d$PValue)
write.table(results.d[o.d,], file="DIclusters.tsv", sep="\t",quote=FALSE, row.names=FALSE)

x=abs(results.d$logFC)
results.d$score = as.integer((1000-0)*((x-min(x))/(max(x)-min(x))) + 0)
results.d[,c("start2", "end2","start1", "end1")] = results.d[,c("start2", "end2","start1", "end1")] + region_start
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
write.table(a, file = "DIclusters2.tsv", sep = "\t",quote = FALSE, 
            row.names = FALSE, col.names = TRUE)
arc=results.d[,c("seqnames2","start2", "end2", "seqnames1","start1", "end1", "logFC")]
write.table(arc, file = paste(region, "di_cluster.arc", sep = ""), sep = "\t",quote = FALSE, 
            row.names = FALSE, col.names = FALSE)
