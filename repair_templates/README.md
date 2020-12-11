# Preparing to generate repair templates

Add nucleotide sequences upstream and downstream of the original genome location used for chimeragenesis.  i.e. when chimerizing PCA1 with CAD2 allele, it will be done around the original PCA1 site. This will require two files, `PCA1.upstream` and `PCA1.downstream`. 

I took the `Genomic DNA +/- 1kb` fasta from [SGD](https://www.yeastgenome.org/locus/S000000499#sequence). The first 1000 nucleotides made `PCA1.upstream` and last 1000 nucleotides made `PCA1.downstream`.

# Copy required files to this directory
```
cp ../codon_shuffle/CAD2.min.fasta ./CAD2.min.cds.fasta
cp ../codon_shuffle/CAD2.low.fasta ./CAD2.low.cds.fasta
cp ../codon_shuffle/CAD2.medium.fasta ./CAD2.medium.cds.fasta
cp ../codon_shuffle/CAD2.high.fasta ./CAD2.high.cds.fasta
cp ../codon_shuffle/CAD2.cds.fasta .
```

## Generate Repair Templates that show the method works
Using longest repair templates (80 bp for each, or 160 bp total) and most diverged CAD2 sequence
```
mkdir -p 01_test_method
python3 ./chimera.py  PCA1 CAD2.min --flanking PCA1 --repair-template-length 160 --unique protein > 01_test_method/PCA1_CAD2.min.RT-160.fasta
python3 ./chimera.py  CAD2.min PCA1 --flanking PCA1 --repair-template-length 160 --unique protein > 01_test_method/CAD2.min_PCA1.RT-160.fasta

# 27 rts each * 2 = 54
```

## Look at variety of repair template lengths and sequence homology
This will be a total of 10 transformations (PCA1-CAD2 and CAD2-PCA1 orientations, with 5 levels of sequence homology)
```
mkdir -p 02_RT_length
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PCA1_CAD2.min.RT-all.fasta  ::: PCA1 ::: CAD2.min :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PCA1_CAD2.low.RT-all.fasta  ::: PCA1 ::: CAD2.low :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PCA1_CAD2.medium.RT-all.fasta  ::: PCA1 ::: CAD2.medium :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PCA1_CAD2.high.RT-all.fasta  ::: PCA1 ::: CAD2.high :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PCA1_CAD2.orig.RT-all.fasta ::: PCA1 ::: CAD2 ::: 40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/CAD2.min_PCA1.RT-all.fasta  ::: CAD2.min ::: PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/CAD2.low_PCA1.RT-all.fasta  ::: CAD2.low ::: PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/CAD2.medium_PCA1.RT-all.fasta  ::: CAD2.medium ::: PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/CAD2.high_PCA1.RT-all.fasta  ::: CAD2.high ::: PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/CAD2.orig_PCA1.RT-all.fasta ::: CAD2 ::: PCA1 ::: 40 44 50 58 68 80 94 110 128 148 160

# 297 RTs each * 10 = 2970
```

## Look at the affects of chimerizing at synonymous sites
```
mkdir -p 03_synonymous_RT
python3 ./chimera.py  PCA1 CAD2.min --flanking PCA1 --repair-template-length 160 --unique dna > 03_synonymous_RT/PCA1_CAD2.min.RT-160-syn.fasta
python3 ./chimera.py  CAD2.min PCA1 --flanking PCA1 --repair-template-length 160 --unique dna > 03_synonymous_RT/CAD2.min_PCA1.RT-160-syn.fasta

# 1217 RTs each * 2 = 2434
```