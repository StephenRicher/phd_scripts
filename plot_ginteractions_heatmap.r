library(ggplot2)
library(ggthemes)
library(plotly)
library(reshape2)
data = read.table('~/x_db/DBuck/s_richer/stephen_test/projects/hic_analysis/post_alignment/hicexplorer/all_regions/GNG12_AS1_DIRAS3/5000/HB2_WT_1-GNG12_AS1_DIRAS3-5000-norm_iced_obs_exp.ginteractions.tsv')


a = dcast(data, V2 ~ V5, value.var = "V7")


p <- ggplot(data, aes(V2, V5)) + geom_tile(aes(fill = V7)) + scale_fill_gradient2(limits=c(-0,2), midpoint = 1)
# scale_fill_gradient2_tableau(palette = "Orange-Blue Diverging", na.value = "grey50", guide = "colourbar")
ggplotly(p)
