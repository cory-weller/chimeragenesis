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
Using longest repair templates (80 bp for each, or 160 bp total) and most diverged CAD2 sequence
```
mkdir -p 01_test_method
python3 ./chimera.py  PCA1 CAD2.55 --flanking PCA1 --repair-template-length 160 --unique protein > 01_test_method/PCA1_CAD2.55.RT-160.fasta
python3 ./chimera.py  CAD2.55 PCA1 --flanking PCA1 --repair-template-length 160 --unique protein > 01_test_method/CAD2.55_PCA1.RT-160.fasta

# 27 rts each * 2 = 54
```

## Look at variety of repair template lengths and sequence homology
This will be a total of 10 transformations (PCA1-CAD2 and CAD2-PCA1 orientations, with 5 levels of sequence homology)
```
mkdir -p 02_RT_length
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein > 02_RT_length/PCA1_CAD2.55.RT-all.fasta  ::: PCA1 ::: CAD2.55 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein > 02_RT_length/PCA1_CAD2.65.RT-all.fasta  ::: PCA1 ::: CAD2.65 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein > 02_RT_length/PCA1_CAD2.74.RT-all.fasta  ::: PCA1 ::: CAD2.74 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein > 02_RT_length/PCA1_CAD2.87.RT-all.fasta  ::: PCA1 ::: CAD2.87 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein > 02_RT_length/PCA1_CAD2.986.RT-all.fasta ::: PCA1 ::: CAD2.986 ::: 40 44 50 58 68 80 94 110 128 148 160

parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein > 02_RT_length/CAD2.55_PCA1.RT-all.fasta  ::: CAD2.55 ::: PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein > 02_RT_length/CAD2.65_PCA1.RT-all.fasta  ::: CAD2.65 ::: PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein > 02_RT_length/CAD2.74_PCA1.RT-all.fasta  ::: CAD2.74 ::: PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein > 02_RT_length/CAD2.87_PCA1.RT-all.fasta  ::: CAD2.87 ::: PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein > 02_RT_length/CAD2.986_PCA1.RT-all.fasta ::: CAD2.986 ::: PCA1 ::: 40 44 50 58 68 80 94 110 128 148 160

# 297 RTs each * 10 = 2970
```

## Look at the affects of chimerizing at synonymous sites
```
mkdir -p 03_synonymous_RT
python3 ./chimera.py  PCA1 CAD2.55 --flanking PCA1 --repair-template-length 160 --unique dna > 03_synonymous_RT/PCA1_CAD2.55.RT-160-syn.fasta
python3 ./chimera.py  CAD2.55 PCA1 --flanking PCA1 --repair-template-length 160 --unique dna > 03_synonymous_RT/CAD2.55_PCA1.RT-160-syn.fasta

# 1193 RTs each * 2 = 2386
```