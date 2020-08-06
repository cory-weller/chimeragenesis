#!/usr/bin/env bash

cat MET12.pep.fasta MET13.pep.fasta | clustalo -i - --outfmt clu -o MET12_MET13.pep.aln

```python
genomic dna = txt[1000:-1000]
```

# cat <(./permute_alignment.py MET12.pep.fasta) MET13.pep.fasta | ./clustalo -i - --outfmt clu
