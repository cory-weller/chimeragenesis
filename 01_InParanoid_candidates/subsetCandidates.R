#!/usr/bin/env Rscript

if(!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}

if(!require(data.table)) {
  install.packages("data.table")
  library(data.table)
}

if(!require(foreach)) {
  install.packages("foreach")
  library(foreach)
}

loadCandidates <- function(candidates_filename) {
  if(! file.exists(candidates_filename)) {
    # load InParanoid table
    inparanoid <- fread(
                        file='sqltable.H.sapiens-S.cerevisiae',
                        fill=TRUE,
                        col.names=c("cluster", "bitscore", "species","inparanoid_score","pdb_id","pidentity")
                      )
    # subset clusters where there is only one S.cerevisiae gene
    dat <- inparanoid[cluster %in% inparanoid[species=="S.cerevisiae"][, .N, by=cluster][N==1, cluster]]
    # leaves 1758 clusters

    # subset clusters where there is only one H.sapiens gene
    dat <- dat[cluster %in% dat[species=='H.sapiens'][, .N, by=cluster][N==1, cluster]]
    # leaves 1266 gene pairs

    # cast from long to wide
    candidates <- dcast(dat, cluster+bitscore~species, value.var="pdb_id")
    fwrite(candidates, file=candidates_filename, quote=F, row.names=F, col.names=T, sep="\t")
  } else {
    cat(paste("loading ", candidates_filename, "\n", sep=""))
    candidates<- fread(candidates_filename)
  }
  return(candidates)
}

loadPhenotypes <- function(phenotypes_filename) {
  if(! file.exists(phenotypes_filename)) {
    auxo <- fread('annotations/auxotrophy_annotations.txt', skip="Gene")
    carbo <- fread('annotations/utilization_of_carbon_source_absent_annotations.txt', skip="Gene")
    respir <- fread('annotations/respiratory_growth_absent_annotations.txt', skip="Gene")
    phenotypes <- rbindlist(list(auxo, carbo, respir))
    fwrite(phenotypes, file=phenotypes_filename, quote=F, row.names=F, col.names=T, sep="\t")
  } else {
    cat(paste("loading ", phenotypes_filename, "\n", sep=""))
    phenotypes <- fread(phenotypes_filename)
  }
  phenotypes[, "Details" := NULL]
  phenotypes[, "Reference" := NULL]
  phenotypes[, "Experiment Type" := NULL]
  phenotypes[, "Experiment Type Category" := NULL]
  return(phenotypes)
}

loadViable <- function(viable_filename) {
  cat(paste("loading ", viable_filename, "\n", sep=""))
  viable <- fread(viable_filename, skip="Gene")
  return(viable)
}

loadStructures <- function(structures_filename) {
  structures <- fread(structures_filename,
                      col.names=c("cluster", "pdb_id", "species", "structure_id",
                                  "method", "resolution", "region"))
  return(structures)
}

setwd("/data/wellerca/yeast-chimeragenesis/01_InParanoid_candidates/")

candidates <- loadCandidates("candidate_pairs.txt")
phenotypes <- loadPhenotypes("phenotypes.txt")

# Fix "BCS1" phenotype label
phenotypes[Gene=="BCS1", Phenotype := "respiratory growth: absent"]

phenotypes.wide <- dcast(phenotypes[, .N, by=list(Gene, Phenotype)], Gene~Phenotype, value.var="N")
phenotypes.wide[, "auxotrophic" := ifelse(is.na(`auxotrophy`), FALSE, TRUE)]
phenotypes.wide[, "non_respiratory" := ifelse(is.na(`respiratory growth: absent`), FALSE, TRUE)]
phenotypes.wide[, "non_carbon_using" := ifelse(is.na(`utilization of carbon source: absent`), FALSE, TRUE)]
phenotypes.wide[, "auxotrophic" := ifelse(is.na(auxotrophy), FALSE, TRUE)]

phenotypes.wide[, c("auxotrophy", "respiratory growth: absent", "utilization of carbon source: absent") := NULL]

# Add Gene IDs to Protein ID table
protein_to_gene_table <- fread("protein_to_gene.txt", fill=TRUE, col.names=c("ProteinID","GeneID"), na.strings=c(""))
setkey(protein_to_gene_table, ProteinID)
candidates[, "H.sapiens_GeneID" := protein_to_gene_table[candidates[, "H.sapiens"], GeneID]]
candidates[, "S.cerevisiae_GeneID" := protein_to_gene_table[candidates[, "S.cerevisiae"], GeneID]]

# remove genes with NAs for either species gene IDs
candidates <- candidates[! is.na(S.cerevisiae_GeneID) & ! is.na(H.sapiens_GeneID)]
# 1213 gene pairs with both Gene IDs

# subset candidates with known auxotrophy / no carbon utilization / no respiration
candidates <- candidates[S.cerevisiae_GeneID %in% phenotypes.wide[,Gene]]
# 296 with above phenotypes

# subset candidates known viable with null
viable <- loadViable("annotations/viable_annotations.txt")
candidates <- candidates[S.cerevisiae_GeneID %in% viable[`Mutant Information`=="null",Gene]]

# 252 viable null candidate gene pairs
structures <- loadStructures('InParanoidStructures.txt')

candidates[, H.sapiens_structure := candidates[,H.sapiens %in% structures[,pdb_id]]]
candidates[, S.cerevisiae_structure := candidates[,S.cerevisiae %in% structures[,pdb_id]]]
candidates[, any_structure := ifelse(S.cerevisiae_structure == TRUE | H.sapiens_structure == TRUE, TRUE, FALSE)]

# Subset candidates with either human OR yeast structure
candidates <- candidates[any_structure==TRUE]

setnames(candidates, "H.sapiens", "H.sapiens_ProteinID")
setnames(candidates, "S.cerevisiae", "S.cerevisiae_ProteinID")

# get auxo phenotypes in table
auxotrophs <- merge(candidates,phenotypes[Phenotype=="auxotrophy"], by.x="S.cerevisiae_GeneID",  by.y="Gene")
fwrite(auxotrophs, file="auxotrophic_candidates.txt",  quote=F, row.names=F, col.names=T, sep="\t")

# merge in phenotypes
candidates <- merge(candidates, phenotypes.wide, by.x="S.cerevisiae_GeneID", by.y="Gene")
fwrite(candidates, file="167-candidates.txt", quote=F, row.names=F, col.names=T, sep="\t")
