#!/usr/bin/env python

import sys
import re
import glob

# filename = sys.argv[1]
file_list = glob.glob("protein_alignments/*.aln")

def unfold_alignment(filename):
    # read text in .aln file
    with open(filename, 'r') as infile:
      text = infile.readlines()[3:]
    
    text = [x.rstrip() for x in text]
    
    # iterate through file and unfold text
    lineCount = 0
    seq1 = ''
    seq2 = ''
    alnScore = ''
    bufferLength = 12
    
    for line in text:
        if lineCount %4 == 0:
            seq1add = line[bufferLength:]
            lineLen = len(seq1add)
            seq1 += seq1add
        elif lineCount %4 == 1:
            seq2add = line[bufferLength:]
            lineLen = len(seq2add)
            seq2 += seq2add
        elif lineCount %4 == 2:
            alnAdd = line[bufferLength:]
            alnScore += alnAdd + ' ' *(lineLen - len(alnAdd))
        elif lineCount %4 == 3:
            pass
        lineCount += 1
    
    # modify score string
    newScore = ''
    for aa1, aa2, score in zip(seq1, seq2,alnScore):
        if aa1 == "-" and aa2 != "-":
            newScore += 's'
        elif aa1 != "-" and aa2 == "-":
            newScore += 'f'
        else:
            newScore += score
    
    return([seq1, seq2, newScore])

###

def get_overall_score(score):
    gaps = re.findall('[fs]+', score)
    nonperfect_align = re.findall('[^\*fs]', score)
    perfect_align = re.findall('\*', score)
    total_gap_length = sum([len(x) for x in gaps])
    total_align_length = len(nonperfect_align) + len(perfect_align)
    overall_score = float(len(perfect_align)) / (len(perfect_align) + len(nonperfect_align) + len(gaps))
    return([overall_score, total_align_length, total_gap_length, len(perfect_align), len(nonperfect_align), len(gaps)])

with open('scores.txt', 'w') as score_file:
    score_file.write('\t'.join(['cluster', 'protein1', 'protein2', 'overall_score', 'total_align_length', 'total_gap_length', 'total_perfect_align', 'total_nonperfect_align', 'n_gaps']) + '\n')
    for filename in file_list:
        cluster, protein1, protein2 = re.split("\.|_", filename.split("/")[1])[0:3]
        with open(filename + ".2", 'w') as outfile:
            seq1, seq2, score = unfold_alignment(filename)
            outfile.write('\n'.join([seq1, seq2, score]))
            outfile.write('\n')
            overall_score, total_align_length, total_gap_length, total_perfect_align, total_nonperfect_align, n_gaps = get_overall_score(score)
            score_file.write('\t'.join([str(x) for x in [cluster, protein1, protein2, overall_score, total_align_length, total_gap_length, total_perfect_align, total_nonperfect_align, n_gaps]]) + '\n')
