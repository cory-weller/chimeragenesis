#!/usr/bin/env Rscript


args <- commandArgs(trailingOnly=TRUE)

fasta_1_filename <- args[1]
fasta_2_filename <- args[2]

fasta_1 <- readLines(fasta_1_filename)
fasta_1 <- fasta_1[2:length(fasta_1)]
fasta_1 <- paste0(fasta_1, collapse="")
fasta_1 <- strsplit(fasta_1, "(?<=.{3})", perl = TRUE)[[1]]


fasta_2 <- readLines(fasta_2_filename)
fasta_2 <- fasta_2[2:length(fasta_2)]
fasta_2 <- paste0(fasta_2, collapse="")
fasta_2 <- strsplit(fasta_2, "(?<=.{3})", perl = TRUE)[[1]]

library(data.table)
library(foreach)

dat <- fread('codons.tab', col.names=c("codon", "AA", "fraction_per_AA", "count_per_thousand", "count"))
dat[, codon := gsub("U", "T", codon)]

# recalculate fraction_per-AA with greater number of digits
dat[, fraction_per_AA := count/sum(count), by=AA]

# Remove codons with fraction_per_AA lower than 10%
setkey(dat, codon, AA)

AAs <- unique(dat[,AA])
codonTable <- unique(dat[, c("codon","AA")])


o <- CJ("codon1"=dat[, codon], "codon2"=dat[, codon])


dat2 <- merge(o, dat, by.x="codon1", by.y="codon")
setnames(dat2, "codon1", "codon")


dat2[, "aa1" := ifelse(substr(codon, 1, 1) == substr(codon2, 1, 1), FALSE, TRUE)]
dat2[, "aa2" := ifelse(substr(codon, 2, 2) == substr(codon2, 2, 2), FALSE, TRUE)]
dat2[, "aa3" := ifelse(substr(codon, 3, 3) == substr(codon2, 3, 3), FALSE, TRUE)]
dat2[, dist := aa1 + aa2 + aa3]
dat2[, c("aa1","aa2","aa3") := NULL]
dat2 <- dat2[dist > 0]

setkey(dat2, codon, AA)
dat <- dat2[fraction_per_AA >= 0.1]


setnames(dat2, "codon", "new_codon")
setnames(dat2, "codon2", "old_codon")



translate <- function(i) {
    return(codonTable[codon==i, AA])
}

# shuffle sequence 2 to be as far away from 1 as possible
codonShuffle <- function(seq_1_codon_static, seq_2_codon_to_shuffle) {
    dt <- dat2[AA == translate(seq_2_codon_to_shuffle) & old_codon == seq_1_codon_static]
    if(nrow(dt) == 0) {
        return(seq_2_codon_to_shuffle)
    } else {
        return(dt[sample(.N, size=1, prob=fraction_per_AA), new_codon])
    }
}

shuffled_seq <- mapply(codonShuffle, seq_1_codon_static = fasta_1, seq_2_codon_to_shuffle = fasta_2)

cat(paste0(shuffled_seq, collapse=""))