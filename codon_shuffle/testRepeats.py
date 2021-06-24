#!/usr/bin/env python

import regex as re
import sys

infileName = sys.argv[1]

with open(infileName, 'r') as infile:
    seq = infile.readlines()[1:]
    seq = ''.join([x.strip() for x in seq])

for k in range(5,50):
    for i in range(len(seq) - k):
        substring = seq[i:(i+k)]
        pattern = re.compile('(%s){s<=1,i<=1,d<=1}' % substring)
        matches = re.findall(pattern, seq, overlapped = False)
        n = len(matches)
        if n > 1:
            print('\t'.join([str(x) for x in [i, k, substring, len(matches)]]))
