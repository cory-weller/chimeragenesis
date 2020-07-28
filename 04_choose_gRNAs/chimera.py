#!/usr/bin/env python

def read_text(filename):
    with open(filename, 'r') as infile:
        rawText = infile.read()
    return(rawText)

def unfold_alignment(rawText):
    text = rawText.splitlines()[3:]
    
    # iterate through file and unfold text
    lineCount = 0
    seq1 = ''
    seq2 = ''
    alnScore = ''
    
    for line in text:
        if lineCount %4 == 0:
            seq1add = line.split()[1]
            lineLen = len(seq1add)
            seq1 += seq1add
        elif lineCount %4 == 1:
            seq2add = line.split()[1]
            lineLen = len(seq2add)
            seq2 += seq2add
        elif lineCount %4 == 2:
            alnAdd = line.lstrip()
            alnPadding = lineLen - len(alnAdd)
            alnScore += ' ' * alnPadding + alnAdd
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


def get_homology_regions(aln_scores, n_perfect_matches):
    boundaries = [(m.start(0), m.end(0)) for m in re.finditer(r"[\*]{%s,}" % n_perfect_matches, aln_scores)]
    return(boundaries)

def get_overall_score(score):
    gaps = re.findall('[fs]+', score)
    nonperfect_align = re.findall(r'[^\*fs]', score)
    perfect_align = re.findall(r'\*', score)
    total_gap_length = sum([len(x) for x in gaps])
    total_align_length = len(nonperfect_align) + len(perfect_align)
    overall_score = float(len(perfect_align)) / (len(perfect_align) + len(nonperfect_align) + len(gaps))
    return([overall_score, total_align_length, total_gap_length, len(perfect_align), len(nonperfect_align), len(gaps)])

# Unnecessary function?
def pad_from_homology(non_homology_regions, nt_padding):
    non_homology_interior = [ (x[0]+nt_padding, x[1]-nt_padding) for x in non_homology_regions[1:-1]]
    non_homology_first = [(non_homology_regions[0][0], non_homology_regions[0][1]-nt_padding)]
    non_homology_last = [(non_homology_regions[-1][0]+nt_padding, non_homology_regions[-1][1])]
    combined_ranges = non_homology_first + non_homology_interior + non_homology_last
    #return(combined_ranges)
    return([x for x in combined_ranges if x[0] <= x[1] and x[0] >= 0])


# with open('scores.txt', 'w') as score_file:
#     score_file.write('\t'.join(['cluster', 'protein1', 'protein2', 'overall_score', 'total_align_length', 'total_gap_length', 'total_perfect_align', 'total_nonperfect_align', 'n_gaps']) + '\n')
#     for filename in file_list:
#         cluster, protein1, protein2 = re.split(r"\.|_", filename.split("/")[1])[0:3]
#         with open(filename + ".2", 'w') as outfile:
#             seq1, seq2, score = unfold_alignment(filename)
#             outfile.write('\n'.join([seq1, seq2, score]))
#             outfile.write('\n')
#             overall_score, total_align_length, total_gap_length, total_perfect_align, total_nonperfect_align, n_gaps = get_overall_score(score)
#             score_file.write('\t'.join([str(x) for x in [cluster, protein1, protein2, overall_score, total_align_length, total_gap_length, total_perfect_align, total_nonperfect_align, n_gaps]]) + '\n')

def open_fasta(filename):
    with open(filename, "r") as infile:
        fasta_text = infile.readlines()
        return(fasta_text)

def get_index_offset(seq):
    seq_offset = []
    index_value = -1
    for i in seq:
        if i != "-":
            index_value += 1
            seq_offset.append(index_value)
        else:
            seq_offset.append("NA")
    return(seq_offset)

def returnCombos(dict1, dict2, ranges):
    for i in ranges:
        for start in range(i[0], i[1]):
            for end in range(start, i[1]):
                yield((d1[start],d2[end]))

def permute_alignment(alignment, n):
    alignment_list = list(alignment)
    for _ in range(n):
        random.shuffle(alignment_list)
        yield(''.join(alignment_list))

def find_matches(string, re_pattern):


if __name__ == "__main__":
    import sys
    import re
    import glob
    import random
    import datatable as dt

    

    print(sys.argv)

    # file_list = glob.glob("protein_alignments/*.aln")
    #args = sys.argv
    args = ["chimera.py", "MET12_MET13.pep.aln", "3", "MET12_flanking.fasta", "MET13_flanking.fasta", "100"]
    clustal_aln_file = args[1]
    n_perfect_matches = int(args[2])
    seq_1_plus_1kb_file = args[3]
    seq_2_plus_1kb_file = args[4]
    permutations = int(args[5])
    patterns = 
    
    clustal_alignment = read_text(clustal_aln_file)
    #print(clustal_alignment)

    unfolded_alignment = unfold_alignment(clustal_alignment)
    #print(unfolded_alignment)

    homology_regions = get_homology_regions(unfolded_alignment[2], n_perfect_matches)
    print(homology_regions)

    non_homology_regions = [(0,homology_regions[0][0])] + [ (homology_regions[x][1], homology_regions[x+1][0]) for x in range(0,len(homology_regions)-1)] + [(homology_regions[-1][1], len(unfolded_alignment[2]))]
    #print(non_homology_regions)


    # Iterate through non-homology-regions, and get all combinations of numbers within



    seq_1_plus_1kb_file = 'MET12_flanking.fasta'
    seq_2_plus_1kb_file = 'MET13_flanking.fasta'

    seq1 = ''.join([x.strip() for x in open_fasta(seq_1_plus_1kb_file)][1:])
    seq2 = ''.join([x.strip() for x in open_fasta(seq_2_plus_1kb_file)][1:])

    # mylist = []
    # for non_homology_range in non_homology_regions:
    #     for i in range(non_homology_range[0], non_homology_range[1]):
    #         seq1_index_i = i - unfolded_alignment[0][:i].count("-")
    #         seq2_index_i = i - unfolded_alignment[1][:i].count("-")
    #         for j in range(i, non_homology_range[1]):
    #             seq1_index_j = j - unfolded_alignment[0][:i].count("-")
    #             seq2_index_j = j - unfolded_alignment[1][:i].count("-")
    #             seq1 [0:i]
    #             seq2 [j:]
    #             mylist.append((i,j))
    # print(len(mylist))



    d1 = get_index_offset(unfolded_alignment[0])
    d2 = get_index_offset(unfolded_alignment[1])
    #c = (  (d1[x],d2[x]) for y in non_homology_regions for x in range(y[0],y[1]))

    c = returnCombos(d1, d2, non_homology_regions)
    print(len(list(c)))

    itertools.permutations(unfolded_alignment[2])

##
a = permute_alignment(unfolded_alignment[2], permutations)
for i in a:
    print(i)
