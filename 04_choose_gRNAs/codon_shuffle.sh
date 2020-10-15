#!/usr/bin/env bash

python3 - << EOF > codons.tab

# Replace table below with codon table from 
# Codon Usage Database: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4932&aa=1&style=N
inputText = """
UUU F 0.59 26.1 (170666)  UCU S 0.26 23.5 (153557)  UAU Y 0.56 18.8 (122728)  UGU C 0.63  8.1 ( 52903)
UUC F 0.41 18.4 (120510)  UCC S 0.16 14.2 ( 92923)  UAC Y 0.44 14.8 ( 96596)  UGC C 0.37  4.8 ( 31095)
UUA L 0.28 26.2 (170884)  UCA S 0.21 18.7 (122028)  UAA * 0.47  1.1 (  6913)  UGA * 0.30  0.7 (  4447)
UUG L 0.29 27.2 (177573)  UCG S 0.10  8.6 ( 55951)  UAG * 0.23  0.5 (  3312)  UGG W 1.00 10.4 ( 67789)

CUU L 0.13 12.3 ( 80076)  CCU P 0.31 13.5 ( 88263)  CAU H 0.64 13.6 ( 89007)  CGU R 0.14  6.4 ( 41791)
CUC L 0.06  5.4 ( 35545)  CCC P 0.15  6.8 ( 44309)  CAC H 0.36  7.8 ( 50785)  CGC R 0.06  2.6 ( 16993)
CUA L 0.14 13.4 ( 87619)  CCA P 0.42 18.3 (119641)  CAA Q 0.69 27.3 (178251)  CGA R 0.07  3.0 ( 19562)
CUG L 0.11 10.5 ( 68494)  CCG P 0.12  5.3 ( 34597)  CAG Q 0.31 12.1 ( 79121)  CGG R 0.04  1.7 ( 11351)

AUU I 0.46 30.1 (196893)  ACU T 0.35 20.3 (132522)  AAU N 0.59 35.7 (233124)  AGU S 0.16 14.2 ( 92466)
AUC I 0.26 17.2 (112176)  ACC T 0.22 12.7 ( 83207)  AAC N 0.41 24.8 (162199)  AGC S 0.11  9.8 ( 63726)
AUA I 0.27 17.8 (116254)  ACA T 0.30 17.8 (116084)  AAA K 0.58 41.9 (273618)  AGA R 0.48 21.3 (139081)
AUG M 1.00 20.9 (136805)  ACG T 0.14  8.0 ( 52045)  AAG K 0.42 30.8 (201361)  AGG R 0.21  9.2 ( 60289)

GUU V 0.39 22.1 (144243)  GCU A 0.38 21.2 (138358)  GAU D 0.65 37.6 (245641)  GGU G 0.47 23.9 (156109)
GUC V 0.21 11.8 ( 76947)  GCC A 0.22 12.6 ( 82357)  GAC D 0.35 20.2 (132048)  GGC G 0.19  9.8 ( 63903)
GUA V 0.21 11.8 ( 76927)  GCA A 0.29 16.2 (105910)  GAA E 0.70 45.6 (297944)  GGA G 0.22 10.9 ( 71216)
GUG V 0.19 10.8 ( 70337)  GCG A 0.11  6.2 ( 40358)  GAG E 0.30 19.2 (125717)  GGG G 0.12  6.0 ( 39359)
"""

inputText = inputText.replace("(","").replace(")","").split()

entries = [inputText[x:(x+5)] for x in range(0, len(inputText), 5)]

for row in entries:
    print('\t'.join([str(x) for x in row]))

EOF

Rscript - PCA1.cds.fasta CAD2.cds.fasta codons.tab << EOF

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

EOF
