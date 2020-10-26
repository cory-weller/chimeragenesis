

`<GENE1>.upstream.fasta` containing nucleotide sequence upstream (5' from perspective of coding strand) of first gene in pair
`<GENE2>.downstream.fasta` containing nucleotide sequence downstream (3' from perspective of coding strand) of second gene in pair


# Codon shuffling CAD2

First, generate codon table
( cd codon_shuffle && if [ ! -f "codons.tab" ]; then python3 printCodons.py > "codons.tab"; fi )

Next, shuffle codons 