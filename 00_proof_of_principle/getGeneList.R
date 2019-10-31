#!/usr/bin/env R
if(! require(rentrez)) {
install.packages("rentrez")
  library(rentrez)
}

if(! require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}

if(! require(data.table)) {
  install.packages("data.table")
  library(data.table)
}

if(! require(foreach)) {
  install.packages("foreach")
  library(foreach)
}

if(! require(doMC)) {
  install.packages("doMC")
  library(doMC)
}

structure_search <- entrez_search(db="structure",
                                  term="Saccharomyces cerevisiae S288C[ORGN]",
                                  retmax=9999,
                                  use_history=FALSE)


o <- foreach(structure_id=structure_search$ids) %do% {
  res <- entrez_summary(db="structure", id=structure_id, rettype="native")
  data.table(t(unlist(res)))
}

saveRDS(o, file="structure_dat.RDS")
# o <- readRDS("structure_dat.RDS")
dat <- rbindlist(o, fill=TRUE)

# get rid of extra organism list columns
extraOrganismListCols <- colnames(dat)[colnames(dat) %like% "organismlist[^$]"]
dat <- dat[, ! ..extraOrganismListCols]
dat[, organismlist := "Saccharomyces cerevisiae S288C"]

# collect and melt pdb accension synonyms
pdb_cols <- colnames(dat)[colnames(dat) %like% "pdbaccsynlist"]
pdb <- dat[, c("uid", "pdbacc", pdb_cols), with=F]
pdb.long <- melt(pdb, measure.vars=pdb_cols)[!is.na(value)][order(uid, pdbacc)][, ! "variable", with=FALSE]
dat[, (pdb_cols) := NULL]

# write pdb ids to file
fwrite(dat[,.(pdbacc)], file="pdb_ids.txt", quote=F, row.names=F, col.names=F, sep="\t")

# convert pdb ids to gid
function get_gid {
  pdbid=${1}
  gid=$(esearch -db protein -query $pdbid | esummary | grep "DocumentSummary><Id" | cut -d ">" -f 3 | cut -d "<" -f 1)
  echo $pdbid $gid
}

export -f get_gid

sed 's/^/get_gid /g' pdb_ids.txt | sed 's/$/ \&\& sleep 0.34/g' > get_gids.sh

nohup bash get_gids.sh > gids.txt &

# now perform hologene lookup using gids




# find homologous genes using HomoloGene

# get homologous human genes
# > entrez_db_searchable("homologene")
# Searchable fields for database 'homologene'
#   ALL    All terms from all searchable fields
#   UID    Unique number assigned to publication
#   FILT   Limits the records
#   TITL   Words in title of publication
#   WORD   Free text associated with record
#   PROP   Properties (formerly Keyword)
#   ORGN   scientific and common names of organism
#   GNID   Gene ID
#   GENE   Gene Name
#   GDSC   Description of gene
#   PUID   protein uids
#   PRAC   protein accessions
#   NUID   nucleotide uids of sequences
#   NCAC   nucleotide accessions of seqeunces
#   UGID   UniGene ID
#   ANCS   scientific and common names of ancestor organism
#   DOM    Domain Name
tpo[gene name] AND human[orgn].
