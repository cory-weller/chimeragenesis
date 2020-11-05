# Preparing to generate repair templates

Add nucleotide sequences upstream and downstream of the original genome location used for chimeragenesis.  i.e. when chimerizing PCA1 with CAD2 allele, it will be done around the original PCA1 site. This will require two files, `PCA1.upstream` and `PCA1.downstream`. 

I took the `Genomic DNA +/- 1kb` fasta from [SGD](https://www.yeastgenome.org/locus/S000000499#sequence). The first 1000 nucleotides made `PCA1.upstream` and last 1000 nucleotides made `PCA1.downstream`.

I created five fasta files for the multiple versions of the CAD2 allele being used for chimeragenesis, Reflecting the % identity shared with `PCA1.cds.fasta` (the original CAD2 allele is 98.6% similar to PCA1).
 * `CAD2.55.cds.fasta`
 * `CAD2.65.cds.fasta`
 * `CAD2.74.cds.fasta`
 * `CAD2.87.cds.fasta`
 * `CAD2.986.cds.fasta`

# Generating repair templates

```

parallel -j 1 python3 ./chimera.py  {1} {2} --flanking PCA1 --repair-template-length {3} ::: PCA1 ::: CAD2 ::: 

```

## Show the method works
Using longest repair templates (80 bp for each, or 160 bp total), both PCA1-CAD2 and CAD2-PCA1 chimeras
```
python3 ./chimera.py  PCA1 CAD2.986 --flanking PCA1 --repair-template-length 160 --unique protein > PCA1_CAD2.986.RT-160.fasta
python3 ./chimera.py  CAD2.986 PCA1 --flanking PCA1 --repair-template-length 160 --unique protein > CAD2.986_PCA1.RT-160.fasta
```

## Look at variety of repair template lengths and sequence homology
```
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein > PCA1_CAD2.55.RT-all.fasta  ::: PCA1 ::: CAD2.55 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein > PCA1_CAD2.65.RT-all.fasta  ::: PCA1 ::: CAD2.65 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein > PCA1_CAD2.74.RT-all.fasta  ::: PCA1 ::: CAD2.74 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein > PCA1_CAD2.87.RT-all.fasta  ::: PCA1 ::: CAD2.87 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein > PCA1_CAD2.986.RT-all.fasta ::: PCA1 ::: CAD2.986 ::: 40 44 50 58 68 80 94 110 128 148 160

parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein > CAD2.55_PCA1.RT-all.fasta  ::: CAD2.55 ::: PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein > CAD2.65_PCA1.RT-all.fasta  ::: CAD2.65 ::: PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein > CAD2.74_PCA1.RT-all.fasta  ::: CAD2.74 ::: PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein > CAD2.87_PCA1.RT-all.fasta  ::: CAD2.87 ::: PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein > CAD2.986_PCA1.RT-all.fasta ::: CAD2.986 ::: PCA1 ::: 40 44 50 58 68 80 94 110 128 148 160
```