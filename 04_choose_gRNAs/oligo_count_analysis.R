#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(foreach)
library(ggthemes)
library(doMC)
registerDoMC(cores=4)

files = list.files(pattern="*.permutations")

o <- foreach(filename=files, .combine="rbind") %dopar% {
    splitname <- strsplit(filename, split="\\.")[[1]]
    specificity <- splitname[1]
    match_length <- splitname[2]
    dat <- fread(filename)
    dat[, "specificity" := specificity]
    dat[, "match_length" := match_length]
    return(dat)
}

o[, specificity := factor(specificity, levels=c("low", "medium", "high"))]

dat.long <- melt(o, measure.vars=c("homology_regions","n_chimeras"))

g <- ggplot(dat.long, aes(x=match_length, y=value, color=specificity)) +
geom_boxplot() +
facet_grid(variable~., scales="free") +
labs(x="Match length", y="Count", title="Permuted MET12/MET13 alignments") +
theme_few(10)

ggsave(g, file="figures/Permuted_alignments.png", width=20, height=20, units="cm")