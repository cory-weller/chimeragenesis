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

##Generate shuffled sequences
```
Rscript shuffle.R PCA1.cds.fasta CAD2.cds.fasta
```
The above script generates four new `fasta` files with various % identity shared with `PCA1.cds.fasta` (in addition to the original `CAD2.cds.fasta`):
 * `CAD2.cds.min.fasta`
 * `CAD2.cds.low.fasta`
 * `CAD2.cds.med.fasta`
 * `CAD2.cds.high.fasta`

 ## Removing homopolymers

 Then, I used `homopolymers.py` to remove homopolymers. See `class fasta` within the script for exact replacements.

 ```
./homopolymers.py CAD2.cds.min.fasta PCA1.cds.fasta > CAD2.min.fasta && rm CAD2.cds.min.fasta
./homopolymers.py CAD2.cds.low.fasta PCA1.cds.fasta > CAD2.low.fasta && rm CAD2.cds.low.fasta
./homopolymers.py CAD2.cds.medium.fasta PCA1.cds.fasta > CAD2.medium.fasta && rm CAD2.cds.medium.fasta
./homopolymers.py CAD2.cds.high.fasta PCA1.cds.fasta > CAD2.high.fasta && rm CAD2.cds.high.fasta

Codon # 278, TTT, is identical to previous codon # 277

$ ./homopolymers.py CAD2.65.fasta
Codon # 219, AAA, is identical to previous codon # 218
Codon # 278, TTT, is identical to previous codon # 277
Codon # 342, AAA, is identical to previous codon # 341
Codon # 343, AAA, is identical to previous codon # 342
Codon # 756, AAA, is identical to previous codon # 755

$ ./homopolymers.py CAD2.74.fasta
Codon # 343, AAA, is identical to previous codon # 342

$ ./homopolymers.py CAD2.87.fasta
Codon # 343, AAA, is identical to previous codon # 342
Codon # 430, AAA, is identical to previous codon # 429
Codon # 659, TTT, is identical to previous codon # 658

$ ./homopolymers.py CAD2.986.fasta
Codon # 659, TTT, is identical to previous codon # 658

 ```

 These stretches were manually edited with the following changes:

 `AAA` -> `AAG`
 `TTT` -> `TTC`
 `GGG` -> `GGT`
 `CCC` -> `CCA`