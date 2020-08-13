#!/usr/bin/env bash

# generate permuted alignment output
permutations="50000"
for length in $(seq 1 5); do
    for specificity in "low" "medium" "high"; do
        python3 chimera.py -p ${permutations} -s ${specificity} -L ${length}  MET12.pep.fasta MET13.pep.fasta
    done
done

# generate true alignment output
for length in $(seq 1 5); do
    for specificity in "low" "medium" "high"; do
        python3 chimera.py -s ${specificity} -L ${length} MET12.pep.fasta MET13.pep.fasta
    done
done

