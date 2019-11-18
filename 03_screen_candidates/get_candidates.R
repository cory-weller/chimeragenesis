#!/usr/bin/env Rscript

library(data.table)
if(!file.exists('Homozygous_diploid_YSC1056_R86-R97.tab')) {
  dat <- fread('Homozygous_diploid_YSC1056_R86-R97.csv', select=1:8)
  fwrite(dat, file='Homozygous_diploid_YSC1056_R86-R97.tab', quote=F, row.names=F, col.names=T, sep="\t")
  rm(dat)
}

KOs <- fread('Homozygous_diploid_YSC1056_R86-R97.tab')
candidates <- fread('../02_score_alignments/167-candidates-scores-names.txt')

inositol_candidates <- fread('../02_score_alignments/51-inositol-candidates-scores-names.txt')
inositol_candidates_KOs <- merge(KOs, inositol_candidates, by.x="ORF", by.y="secondaryIdentifier")
inositol_candidates_KOs <- inositol_candidates_KOs[Comment==""][Plate != 71 & Plate != 72][!is.na(`Record number`)][Batch != ""][!duplicated(ORF)]

# left with 41 inositol auxotrophs

# get non-respiratory (can only grow on YPD, not YPGlycerol)
non_respiratory_candidates_KOs <- merge(KOs, candidates[non_respiratory==TRUE & non_carbon_using==FALSE & auxotrophic==FALSE], by.x="ORF", by.y="secondaryIdentifier")
non_respiratory_candidates_KOs <- non_respiratory_candidates_KOs[Comment==""][Plate != 71 & Plate != 72][!is.na(`Record number`)][Batch != ""][!duplicated(ORF)]

# left with 21 non-respiratory 
merged_candidates <- merge(KOs, candidates, by.x="ORF", by.y="secondaryIdentifier")[,c("ORF","name")]
fwrite(merged_candidates, file='candidates_with_KO_locations', quote=F, row.names=F, col.names=T, sep="\t")


table1 <- non_respiratory_candidates_KOs[,c("ORF","Plate","Row","Col","overall_score","name")][order(Plate)]
table1[,type := "non_respiratory"]


table2 <- inositol_candidates_KOs[,c("ORF","Plate","Row","Col","overall_score","name")][order(Plate)]
table2[,type := "inositol_auxotroph"]

combined <- rbindlist(list(table1,table2))

setkey(combined, Plate, Row, Col)

fwrite(combined, file="63_phenotype_candidates.txt", quote=F, row.names=F, col.names=T, sep="\t")
