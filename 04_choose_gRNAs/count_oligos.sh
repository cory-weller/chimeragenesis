#!/usr/bin/env bash

permutations="50000"

for length in $(seq 1 5); do
    for specificity in "low" "medium" "high"; do
        python3 chimera.py -p ${permutations} -s ${specificity} -L ${length}  MET12.pep.fasta MET13.pep.fasta
    done
done

