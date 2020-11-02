# README

## Download PCA1 and CAD2 coding sequence

PCA1 retrieved from https://www.yeastgenome.org/locus/S000000499

Renamed from `S288C_YBR295W_PCA1_coding.fsa` to `PCA1.cds.fasta`

CAD2 allele retrieved from https://www.ebi.ac.uk/ena/browser/view/AB027571

Renamed from `AB027571.1.fasta` to `CAD2.cds.fasta`

Codon table retrieved from https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4932&aa=1&style=N

Body of table's data saved as text file `codonTable.txt` within `


```
mv S288C_YBR295W_PCA1_coding.fsa PCA1.cds.fasta
mv AB027571.1.fasta CAD2.cds.fasta
```

( cd ../codon_shuffle && if [ ! -f "codons.tab" ]; then python3 printCodons.py > "codons.tab"; fi )
