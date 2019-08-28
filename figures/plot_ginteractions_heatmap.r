library(ggplot2)
library(ggthemes)
library(plotly)
library(reshape2)
data = read.table('~/x_db/DBuck/s_richer/stephen_test/projects/hic_analysis/post_alignment/hicexplorer/all_regions/GNG12_AS1_DIRAS3/5000/HB2_WT_1-GNG12_AS1_DIRAS3-5000-norm_iced_obs_exp.ginteractions.tsv')
data$V7 = ifelse(data$V7 > 2, 2, data$V7)

data = data[data$V2 > 67000000 & data$V5 <  69000000,]

a = dcast(data, V2 ~ V5, value.var = "V7")
rownames(a) = a[,1]
a = as.matrix(a[,-1])

p <- plot_ly(x = rownames(a), y = colnames(a), z = t(a), colors = "PiYG", type = "heatmap")
p
