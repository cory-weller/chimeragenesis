# README

## Download PCA1 and CAD2 coding sequence

PCA1 retrieved from https://www.yeastgenome.org/locus/S000000499

Renamed from `S288C_YBR295W_PCA1_coding.fsa` to `PCA1.cds.fasta`

CAD2 allele retrieved from https://www.ebi.ac.uk/ena/browser/view/AB027571

Renamed from `AB027571.1.fasta` to `CAD2.cds.fasta`:

```
mv S288C_YBR295W_PCA1_coding.fsa PCA1.cds.fasta
mv AB027571.1.fasta CAD2.cds.fasta
```

Codon table retrieved from https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4932&aa=1&style=N

Body of the above table's data saved as text file `codonTable.txt` 

The contents of `codonTable.txt` is parsed and processed by `printCodons.py` into `codons.tab`:

```
( if [ ! -f "codons.tab" ]; then python3 formatCodons.py > "codons.tab"; fi )
```

Generate shuffled sequences
```
Rscript shuffle.R PCA1.cds.fasta CAD2.cds.fasta
cat *.tmp | fold > CAD2.shuffled.fasta && rm *.tmp
```
