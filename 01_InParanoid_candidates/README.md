# Retrieve InParanoid data

```
wget http://inparanoid.sbc.su.se/download/8.0_current/Orthologs_other_formats/H.sapiens/InParanoid.H.sapiens-S.cerevisiae.tgz
# InParanoid.H.sapiens-S.cerevisiae.tgz md5: a4f07d1e2b1a162700ba8113c9470395

# un-archive
tar -zxvf InParanoid.H.sapiens-S.cerevisiae.tgz
```

# Get Uniprot data for each gene in InParanoid table
```
# run script to pull uniprot structure data
sbatch get_uniprot_data.slurm

# create zip archive from uniprot data
cd ./uniprot_db/ && \
  zip -Tm InParanoidGenes-uniprot.zip *.uniprot && \
  mv InParanoidGenes-uniprot.zip ../ && \
  cd ../ && \
  rmdir uniprot_db

# can now access uniprot data via unzip -p e.g.:
# unzip -p InParanoidGenes-uniprot.zip ${uniprotid}.uniprot
```

# extract structure information for genes in InParanoid table

```
while read cluster bitscore species inparanoidscore pdb_id pidentity; do
  unzip -p InParanoidGenes-uniprot.zip ${pdb_id}.uniprot | \
  grep -E "^DR\s+PDB;" | \
  sed -r 's/DR\s+PDB;\s+//g' | \
  sed 's/; /\t/g' | \
  sed "s/^/${cluster}\t${pdb_id}\t${species}\t/g"
done < sqltable.H.sapiens-S.cerevisiae > InParanoidStructures.txt
```

The result should be a 7-column tab-delimited file `InParanoidStructures.txt`
```
1       Q6P2Q9  H.sapiens       3E9L    X-ray   1.95 A  A=1760-2016.
1       Q6P2Q9  H.sapiens       3ENB    X-ray   1.85 A  A/B=1769-1990.
1       Q6P2Q9  H.sapiens       3JCR    EM      7.00 A  A=1-2335.
1       Q6P2Q9  H.sapiens       3LRU    X-ray   1.85 A  A/B=1831-1990.
1       Q6P2Q9  H.sapiens       4JK7    X-ray   1.40 A  A/B=1769-1990.
1       Q6P2Q9  H.sapiens       4JK8    X-ray   1.15 A  A/B=1769-1990.
```

With columns that correspond to

| x | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
| - | - | - | - | - | - | - | - |
| field | cluster | Protein (uniprot) | Species | Structure (pdb) | Method | Resolution (A) | Region (fragment=start-stop) |
| example | 1 | Q6P2Q9 | H.sapiens | 3E9L | X-ray | 1.95 A | A=1760-2016. |

# create table to convert protein ID to gene ID
```
while read cluster bitscore species inparanoidscore pdb_id pidentity; do
  geneID=$(unzip -p InParanoidGenes-uniprot.zip ${pdb_id}.uniprot | \
  grep -E "GN\\s+Name=" | \
  head -n 1 | \
  cut -d "=" -f 2 | \
  cut -d ";" -f 1 | \
  cut -d " " -f 1)
  echo -e "${pdb_id}\t${geneID}"
done < sqltable.H.sapiens-S.cerevisiae > protein_to_gene.txt
```

# Retrieve Saccharomyces Genome Database phenotypes annotations

Download the following:
  * [`auxotrophy_annotations.txt`](https://www.yeastgenome.org/phenotype/auxotrophy)
  * [`utilization_of_carbon_source_absent_annotations.txt`](https://www.yeastgenome.org/phenotype/absent_utilization_of_carbon_source)
  * [`respiratory_growth_absent_annotations.txt`](https://www.yeastgenome.org/phenotype/absent_respiratory_growth)
  * [`viable_annotations.txt`](https://www.yeastgenome.org/phenotype/viable)



# RAD26
# 126 Q03468    H.sapiens 4CVO X-ray 1.85 A            A=84-160.
# 126 P40352 S.cerevisiae 5VVR    EM 5.80 A            M=1-1085.
