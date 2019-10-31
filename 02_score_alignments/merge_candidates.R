#!/usr/bin/env Rscript

if(!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}

if(!require(data.table)) {
  install.packages("data.table")
  library(data.table)
}

scores <- fread('scores.txt')
candidates <- fread('../01_InParanoid_candidates/167-candidates.txt')
dat <- merge(scores, candidates, by="cluster")

# dat[, sgd_link := paste("https://www.yeastgenome.org/locus/", S.cerevisiae_GeneID, sep="")]

fwrite(dat, file="167-candidates-scores.txt", sep="\t", quote=F, row.names=F, col.names=T)

# in terminal, open list of S.cerevisiae geneIDs with
# cut -d $'\t' -f 10 167-candidates-scores.txt | sed 's/^/"/g' | sed 's/$/"/g' | gedit -
# and paste into SGD https://yeastmine.yeastgenome.org/yeastmine/bag.do
# to create summary, save as summary.csv and upload to this directory

sgd_summary <- fread('summary.csv')

sgd_summary <- sgd_summary[reason=="MATCH", c("symbol","name","length","secondaryIdentifier","primaryIdentifier")]

sgd_summary <- merge(dat, sgd_summary, by.x="S.cerevisiae_GeneID", by.y="symbol")

sgd_summary[, hyperlink := paste("https://www.yeastgenome.org/locus/", primaryIdentifier, sep="")]
setcolorder(sgd_summary, c("cluster","protein1","protein2","overall_score","total_align_length","total_gap_length","total_perfect_align","total_nonperfect_align","n_gaps","bitscore","H.sapiens_ProteinID","S.cerevisiae_ProteinID","H.sapiens_GeneID","S.cerevisiae_GeneID","H.sapiens_structure","S.cerevisiae_structure","name","length","secondaryIdentifier","primaryIdentifier","hyperlink"))
fwrite(sgd_summary, file="167-candidates-scores-names.txt", sep="\t", quote=F, row.names=F, col.names=T)




dat[H.sapiens_structure==TRUE & S.cerevisiae_structure==TRUE]

ggplot(data=dat, mapping=aes(x=overall_score)) + geom_histogram(bins=50)


# good examples with score > 0.5
# P07259
# ----------MRLLGAA-------AV-AALGRG---RAP-ASLGWQRKQ--VNWKACRWSSSGVIPNEKIRNIGISAHIDSGKTTLTERVLYYTGRIAKMHEVKGKDGVGAVMDSMELERQRGITI
# MSVQKMMWVPRKMVGGRIPFFTCSKVFSGFSRRSFHESPLARSTYEEEKVLVDEIKQKLTPDDIGRCNKLRNIGISAHIDSGKTTFTERVLYYTKRIKAIHEVRGRDNVGAKMDSMDLEREKGITI
# ssssssssss :::*. sssssss *s:.:.* sss.:*s*   ::.::ss*:    : : ..:   :*:***************:******** **  :***:*:*.*** ****:***::****

# P19414
# MAPYSLLVTR-LQKALGVRQYHVASVLCQRAKVAMSHFEPNEYIHYDLLEKNINIVRKRLNRPLTLSEKIVYGHLDDPASQEIERGKSYLRLRPDRVAMQDATAQMAMLQFISSGLSKVAVPSTIH
# -----MLSARSAIKRPIVRGLATVSNLTRDSKVNQNLLEDHSFINYKQNVETLDIVRKRLNRPFTYAEKILYGHLDDPHGQDIQRGVSYLKLRPDRVACQDATAQMAILQFMSAGLPQVAKPVTVH
# fffff:* :*s  *   **   ..* * : :**  . :* :.:*:*.   :.::*********:* :***:******* .*:*:** ***:******* ********:***:*:** :** * *:*

# # bad examples with score < 0.3
# P34756
# -------------------------------------------------------MATDDKTSPTLDSANDLPRSPT---------SPSHLTHFKP------------LTPDQDEPPFKSAYSSFV
# MSSEEPHASISFPDGSHVRSSSTGTSSVNTIDATLSRPNYIKKPSLHIMSTSTTSTTTDLVTNPILSNISVPKISPPTSSSIATATSTSHVTGTASHSNIKANANTSTSVNKKNLPPTTSGRIPSS
# sssssssssssssssssssssssssssssssssssssssssssssssssssssss :**  *.* *.. .    ** sssssssss* **:*    ssssssssssss . .:: ** .*.     

# Q06706
# -----------------MRNLKLFRTLEFRDIQ----GPGNPQCFSLR-------TEQ-GTVLIGSEHGLIEVDPVSREVKNEVSLVAEGFLPEDGSGRIVGVQDLLDQESVCVATASGDVILC--
# MVEHDKSGSKRQELRSNMRNLITLNKGKFKPTASTAEGDEDDLSFTLLDSVFDTLSDSITCVLGSTDIGAIEVQQFMKDGSRN---VLASFNIQTFDDKLLSFVHFADINQLVFVFEQGDIITATY
# sssssssssssssssss****  :.. :*:   ssss*  :  .*:* sssssss::.s  ** .:: * ***: . :: ..:fff*  .*  :  ..:::.. .: * :.: ..  .**:* .ss
