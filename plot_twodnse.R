#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)

round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}

localMinima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(.Machine$integer.max, x)) < 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  # Return array of local minima indexes
  #y
  a = logical(length(x))
  a[y[1]] = TRUE
  a
}

# Paper - https://www.nature.com/articles/s41467-018-05691-7

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop("One argument must be supplied (input file).n", call.=FALSE)
}
dir = dirname(args[1])
if (dir == ".") {
  dir = "./"
}

twonsd = read.table(args[1], sep = '\t', col.names = c("sample", "region", "binsize", "twodnSE"))
twonsd$sample = sapply(strsplit(as.character(twonsd$sample), "-"), `[`, 1)

twonsd = twonsd %>%
  group_by(sample, region) %>%
  mutate(minima = localMinima(twodnSE))


mean_bin = round_any(mean(twonsd[twonsd$minima == TRUE, "binsize"][[1]]), 1000)
region = unique(twonsd$region)

g = ggplot(twonsd, aes(x = binsize, y = twodnSE, colour = sample)) +
  geom_line() + geom_point()  +
  geom_label_repel(aes(label = ifelse(minima == TRUE, as.character(binsize), '')), 
                   colour = 'black', min.segment.length = 0) +
  geom_vline(xintercept = mean_bin) +
  theme_pubr(legend = "right") +
  scale_colour_brewer(name = NULL, palette = "Dark2") +
  labs(title = NULL,
       subtitle = paste(region, ": optimal binsize =", mean_bin),
       tag = NULL,
       x = "Bin size",
       y = "2-D normalized structural entropy",
       caption = NULL)
g
width = 10
height = (9/16) * width
ggsave(filename = paste(dir,region, "-2DnSE.png", sep = ""), plot = g, 
       dpi = 300, width = width, height = height)
