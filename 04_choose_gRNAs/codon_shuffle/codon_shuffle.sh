#!/usr/bin/env bash

if [ ! -f "codons.tab" ]; then
    python3 printCodons.py > "codons.tab"
fi


Rscript shuffle.R PCA1.cds.fasta CAD2.cds.fasta codons.tab 