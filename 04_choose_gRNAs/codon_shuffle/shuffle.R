#!/usr/bin/env Rscript


args <- commandArgs(trailingOnly=TRUE)


# If no args provided, treat as interactive session and provide args here
if(length(args) == 1) {
    args <- c("CAD2.cds.fasta", "PCA1.cds.fasta", "codon_shuffle/codons.tab")
}

fasta_1_filename <- args[1]
fasta_2_filename <- args[2]
codons_filename <- args[3]



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

dat <- fread(codons_filename, col.names=c("codon", "AA", "fraction_per_AA", "count_per_thousand", "count"))
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

setkey(dat2, codon, AA)
dat <- dat2[fraction_per_AA >= 0.1]


setnames(dat2, "codon", "new_codon")
setnames(dat2, "codon2", "old_codon")

translate <- function(i) {
    return(codonTable[codon==i, AA])
}

# shuffle sequence 2 to be <method> divergent from 1.
# method = "max" will return maximally divergent sequence
# method = "min" will return minimally divergent sequence (i.e. as similar as possible on nucleotide level)
# method = "random" 
codonShuffle <- function(seq_1_codon_static, seq_2_codon_to_shuffle, method) {
    dt <- dat2[AA == translate(seq_2_codon_to_shuffle) & old_codon == seq_1_codon_static]
    if(nrow(dt) == 0) {
        return(seq_2_codon_to_shuffle)
    } else {
        if(method == "max") {
            return(dt[which(dist==max(dist))][sample(.N, size=1, prob=fraction_per_AA), new_codon])
        } else if(method == "high") {
            return(dt[sample(.N, size=1, prob=((1 + dist)^1.7 * fraction_per_AA)), new_codon])
        } else if(method == "min") {
            return(dt[which(dist==min(dist))][sample(.N, size=1, prob=fraction_per_AA), new_codon])
        } else if(method == "low") {
            return(dt[sample(.N, size=1, prob=((4 - dist)^4 * fraction_per_AA)), new_codon])
        } else if(method == "random") {
            return(dt[sample(.N, size=1, prob=fraction_per_AA), new_codon])
        } else if(method == "equal") {
            return(dt[sample(.N, size=1), new_codon])
        } else {
            stop("Invalid method. Use method = 'max', 'min' or 'random'")
        }
    }
}

getStringDistance <- function(a, b) {
    a <- strsplit(a, split = "")[[1]]
    b <- strsplit(b, split = "")[[1]]
    # Returns length of differences between A and B if they are of same length
    if(length(a) != length(b)) {
        stop("strings are not of equal length!")
    } else {
        return(sum(a != b))
    }
}


shuffled_random <- paste0(mapply(codonShuffle, seq_1_codon_static = fasta_1, seq_2_codon_to_shuffle = fasta_2, method="random"), collapse="")
shuffled_max <- paste0(mapply(codonShuffle, seq_1_codon_static = fasta_1, seq_2_codon_to_shuffle = fasta_2, method="max"), collapse="")
shuffled_min <- paste0(mapply(codonShuffle, seq_1_codon_static = fasta_1, seq_2_codon_to_shuffle = fasta_2, method="min"), collapse="")
shuffled_low <- paste0(mapply(codonShuffle, seq_1_codon_static = fasta_1, seq_2_codon_to_shuffle = fasta_2, method="low"), collapse="")
shuffled_high <- paste0(mapply(codonShuffle, seq_1_codon_static = fasta_1, seq_2_codon_to_shuffle = fasta_2, method="high"), collapse="")
shuffled_equal <- paste0(mapply(codonShuffle, seq_1_codon_static = fasta_1, seq_2_codon_to_shuffle = fasta_2, method="equal"), collapse="")

fasta_1_collapsed <- paste0(fasta_1, collapse="")

n_min <- getStringDistance(shuffled_min, fasta_1_collapsed)
n_low <- getStringDistance(shuffled_low, fasta_1_collapsed)
n_random <- getStringDistance(shuffled_random, fasta_1_collapsed)
n_high <- getStringDistance(shuffled_high, fasta_1_collapsed)
n_max <- getStringDistance(shuffled_max, fasta_1_collapsed)

plot(x=0:4, y=c(n_min, n_low, n_random, n_high, n_max))
abline(a=(1633/4))

getStringDistance(shuffled_equal, fasta_1_collapsed)

cat(paste0(shuffled_seq, collapse=""))

A to (not G) and C to (not T)
shuffled_low <- paste0(mapply(codonShuffle, seq_1_codon_static = fasta_1, seq_2_codon_to_shuffle = fasta_2, method="low"), collapse="")
getStringDistance(shuffled_low, fasta_1_collapsed)
shuffled_high <- paste0(mapply(codonShuffle, seq_1_codon_static = fasta_1, seq_2_codon_to_shuffle = fasta_2, method="high"), collapse="")
getStringDistance(shuffled_high, fasta_1_collapsed)


shuffled_random <- paste0(mapply(codonShuffle, seq_1_codon_static = fasta_1, seq_2_codon_to_shuffle = fasta_2, method="random"), collapse="")
getStringDistance(shuffled_random, fasta_1_collapsed)
