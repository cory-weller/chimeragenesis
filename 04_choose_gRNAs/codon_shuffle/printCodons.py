#!/usr/bin/env python

# reformats the codon table codonTable.txt, from
# the Codon Usage Database: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4932&aa=1&style=N

with open("codonTable.txt", "r") as infile:
    inputText = infile.read().replace("(","").replace(")","").split()

entries = [inputText[x:(x+5)] for x in range(0, len(inputText), 5)]

# Returns the distance and number of bad transitions between two codons
# Specifically, transitioning FROM codon 1 TO codon2
def getDistance(codon1, codon2):
    distance = 0
    badTransitions = 0
    for i in zip(codon1, codon2):
        if i[0] != i[1]:
            distance += 1
        if i[0] == 
    return(distance)


codons = []

for row in entries:
    print('\t'.join([str(x) for x in row]))
    codon, AA = row[0:2]
    codons[codon] = AA



