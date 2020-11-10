# Generating codon-shuffled versions of the CAD2 allele

## Retrieving PCA1 and CAD2 coding sequences

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

The contents of `codonTable.txt` is parsed and processed by [`printCodons.py`](printCodons.py) into `codons.tab`:

```
( if [ ! -f "codons.tab" ]; then python3 formatCodons.py > "codons.tab"; fi )
```

## Generate shuffled sequences
```
Rscript shuffle.R PCA1.cds.fasta CAD2.cds.fasta
```
The script [`shuffle.R`](shuffle.R)generates four new `fasta` files with various % identity shared with `PCA1.cds.fasta` (in addition to the original `CAD2.cds.fasta` which is 3651 nt long, has a Levenshtein distance of 52 for 98.58% identity with PCA1):
 * `CAD2.cds.min.fasta`
 * `CAD2.cds.low.fasta`
 * `CAD2.cds.med.fasta`
 * `CAD2.cds.high.fasta`

 ## Removing homopolymers

 Then, I used [`homopolymers.py`](homopolymers.py) to remove homopolymers. See `class fasta` within the script for exact replacements.

```
./homopolymers.py CAD2.cds.min.fasta PCA1.cds.fasta > CAD2.min.fasta && rm CAD2.cds.min.fasta
./homopolymers.py CAD2.cds.low.fasta PCA1.cds.fasta > CAD2.low.fasta && rm CAD2.cds.low.fasta
./homopolymers.py CAD2.cds.medium.fasta PCA1.cds.fasta > CAD2.medium.fasta && rm CAD2.cds.medium.fasta
./homopolymers.py CAD2.cds.high.fasta PCA1.cds.fasta > CAD2.high.fasta && rm CAD2.cds.high.fasta
```

The following code block will show presence of 5+ repeat homopolymers in the original `CAD2.cds.fasta`:
```
for file in CAD2.cds.fasta; do
    for nucleotide in A C T G; do
        cat ${file} | tr -d "\n" | grep -E "[${nucleotide}]{5,}"
    done
done
```

Whereas when checking the newly generated `CAD2` files, no homopolymers are found:
```
for file in CAD2.high.fasta CAD2.medium.fasta CAD2.low.fasta CAD2.min.fasta; do
    for nucleotide in A C T G; do
        cat ${file} | tr -d "\n" | grep -E "[${nucleotide}]{5,}"
    done
done
```

## Final levels of identity for variable CAD2
| CAD2 version                             | Percent identity| Levenshtein distance |
|------------------------------------------|-----------------|----------------------|
| [`CAD2.cds.fasta`](CAD2.cds.fasta)       | 98.58%          |  52                  |
| [`CAD2.high.fasta`](CAD2.high.fasta)     | 86.88%          |  479                 |
| [`CAD2.medium.fasta`](CAD2.medium.fasta) | 74.23%          |  941                 |
| [`CAD2.low.fasta`](CAD2.low.fasta)       | 65.38%          |  1264                |
| [`CAD2.min.fasta`](CAD2.min.fasta)       | 55.49%          |  1625                |

Alignment of all protein sequences [`CAD2.pep.aln`](CAD2.pep.aln) shows expected 100% identity.
