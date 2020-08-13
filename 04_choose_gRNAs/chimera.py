#!/usr/bin/env python

####################
# DEFINE FUNCTIONS #
####################

def read_text(filename):
    with open(filename, 'r') as infile:
        rawText = infile.read()
    if debug:
        print("Read in text:\n" + rawText)
    return rawText

def unfold_alignment(rawText):
    text = rawText.splitlines()[3:]
    
    lineCount = 0
    seq1 = ''
    seq2 = ''
    alnScore = ''
    
    for line in text:
        if lineCount %4 == 0:
            seq1add = line.split()[1]
            line_len = len(seq1add)
            seq1 += seq1add
        elif lineCount %4 == 1:
            seq2add = line.split()[1]
            line_len = len(seq2add)
            seq2 += seq2add
        elif lineCount %4 == 2:
            alnAdd = line.lstrip()
            alnPadding = line_len - len(alnAdd)
            alnScore += ' ' * alnPadding + alnAdd
        elif lineCount %4 == 3:
            pass
        lineCount += 1
    # modify score string
    new_score = ''
    for aa1, aa2, score in zip(seq1, seq2,alnScore):
        if aa1 == "-" and aa2 != "-":
            new_score += 's'
        elif aa1 != "-" and aa2 == "-":
            new_score += 'f'
        else:
            new_score += score
    return [seq1, seq2, new_score]

def get_homology_regions(alignment_scores, length, specificity):
    if specificity == "high":
        pattern = r"[*]"
    elif specificity == "medium":
        pattern = r"[*;]"
    elif specificity == "low":
        pattern = r"[*;.]"
    pattern += "{%s,}" % length
    boundaries = [(m.start(0), m.end(0)) for m in re.finditer(pattern, alignment_scores)]
    return boundaries

def get_non_homology_regions(homology_regions, alignment_length):
    return ( [(0,homology_regions[0][0])] +
                [ (homology_regions[x][1], homology_regions[x+1][0]) for x in range(0,len(homology_regions)-1)] +
                [(homology_regions[-1][1], alignment_length)] )

def get_overall_score(score):
    gaps = re.findall('[fs]+', score)
    nonperfect_align = re.findall(r'[^\*fs]', score)
    perfect_align = re.findall(r'\*', score)
    total_gap_length = sum([len(x) for x in gaps])
    total_align_length = len(nonperfect_align) + len(perfect_align)
    overall_score = float(len(perfect_align)) / (len(perfect_align) + len(nonperfect_align) + len(gaps))
    return [overall_score, total_align_length, total_gap_length, len(perfect_align), len(nonperfect_align), len(gaps)]

def open_fasta(filename):
    with open(filename, "r") as infile:
        fasta_text = infile.readlines()
        return fasta_text

def get_index_offset(seq):
    seq_offset = []
    index_value = -1
    for i in seq:
        if i != "-":
            index_value += 1
            seq_offset.append(index_value)
        else:
            seq_offset.append("NA")
    return seq_offset

def return_combos(dict1, dict2, ranges):
    for i in ranges:
        for start in range(i[0], i[1]):
            for end in range(start, i[1]):
                value1 = dict1[start]
                value2 = dict2[end]
                if not (value1 == 'NA' or value2 == 'NA'):
                    yield(dict1[start], dict2[end])
                else:
                    continue

def permute_alignment(alignment, n):
    alignment_list = list(alignment)
    for _ in range(n):
        random.shuffle(alignment_list)
        yield(''.join(alignment_list))

def wrap_fasta(seq, line_length=60):
    return '\n'.join([seq[i:i+line_length] for i in range(0, len(seq), line_length)])

def permute_peptide_fasta(fasta_text):
    split_text = fasta_text.rstrip().split('\n')
    header = split_text[0]
    seq = list(''.join(split_text[1:]).rstrip('*'))
    random.shuffle(seq)
    permuted_seq = header + '\n' + wrap_fasta(''.join(seq) + '*') + '\n'
    if debug == True:
        print(permuted_seq)
    return permuted_seq

def run_clustal(b_alignment):
    output = subprocess.run(["./clustalo", "-i", "-", "--outfmt", "clu"], input=b_alignment, stdout=PIPE).stdout.decode()
    if debug:
        print(output)
    return output

def run_permutations(n_permutations, pep1, pep2, length, specificity):
    for _ in range(n_permutations):
        yield permuted_alignment(random.randint(0,1), pep1, pep2, length, specificity)

# Unused
# def write_permutations(permutations, n_permutations, outfile_name="permutations.out"):
#     print("Writing %s permutations to %s" % (n_permutations, outfile_name))
#     with open(outfile_name, "w") as outfile:
#         current_i = 0
#         for i in permutations:
#             current_i += 1
#             if current_i % 10 == 0:
#                 print(current_i)
#             outfile.write('\n'.join(unfold_alignment(i)) + '\n\n')

##################
# DEFINE CLASSES #
##################

class permuted_alignment:
    help = 'clustal alignment of permuted protein sequence'
    def __init__(self, zero_or_one, pep1, pep2, length, specificity):
        if zero_or_one == 0:
            self.pep1 = pep1
            self.pep2 = permute_peptide_fasta(pep2)
        elif zero_or_one == 1:
            self.pep1 = permute_peptide_fasta(pep1)
            self.pep2 = pep2
        self.clustal = run_clustal(
                                    (self.pep1 + self.pep2).encode() 
                                    )
        self.aln1, self.aln2, self.scores = unfold_alignment(self.clustal)  # unfolded clustal alignment
        self.homology_regions = get_homology_regions(self.scores, length, specificity)
        self.non_homology_regions = get_non_homology_regions(self.homology_regions, len(self.aln1))
        self.offset1 = get_index_offset(self.aln1)
        self.offset2 = get_index_offset(self.aln2)
        self.combos = list(return_combos(self.offset1, self.offset2, self.non_homology_regions))

#######################################################################################################################
#                                                        MAIN
#######################################################################################################################

if __name__ == "__main__":

    ##################
    # IMPORT MODULES #
    ##################

    import os
    import sys
    import re
    import random
    import argparse
    import subprocess
    from subprocess import PIPE

    parser = argparse.ArgumentParser()
    parser.add_argument('pep1_filename', metavar='pep1.fasta', type=str,
                        help='First peptide fasta file')
    parser.add_argument('pep2_filename', metavar='pep2.fasta', type=str,
                        help='Second peptide fasta file')
    parser.add_argument("-d", "--debug", 
                        help="Assist with debugging by increasing output verbosity",
                        action="store_true")
    parser.add_argument("-L", "--length",
                        type=int,
                        nargs='?',
                        const=1,
                        default=4,
                        help='''Minimum number of consecutive amino acids above the
                                specificity threshold (defined by -S)''')
    parser.add_argument("-s", "--specificity",
                        type=str,
                        nargs='?',
                        const=1,
                        default="high",
                        help=("Degree of matching in clustal alignment output." +
                              "Used with -L to determine minimum length of homologous region.\n" +
                              "low includes [.:*], medium includes [:*], high includes [*]."))
    parser.add_argument("-p", "--permutations",
                        type=int,
                        nargs='?',
                        const=1,
                        help = ("Number of alignment permutations to run. Default: 100" +
                                "Used with -L to determine minimum length of homologous region.\n" +
                                "low includes [.:*], medium includes [:*], high includes [*]."))
    args = parser.parse_args()
    if args.debug:
        debug = True
        print("debug mode on, increasing verbosity.")
    else:
        debug = False

    print("""running with length = %s and specificity = %s""" % (args.length, args.specificity))

#######################################################################################################################

    ######################
    # VALIDATE ARGUMENTS #
    ######################

    try: assert args.length > 0, "ERROR: length must be > 0"
    except AssertionError as error: sys.exit(error)

    try: assert args.specificity in ["low", "medium", "high"], "specificity must be low, medium, or high"
    except AssertionError as error: sys.exit(error)
    
#######################################################################################################################

    ##################
    # VALIDATE INPUT #
    ##################

    missing_files = 0
    for filename in [args.pep1_filename, args.pep2_filename]:
        if not os.path.isfile(filename):
            print("Error: File %s does not exist" % filename)
            missing_files += 1
    if not os.path.isfile("clustalo"):
        print("Error: clustalo binary not found")
        missing_files += 1
    if missing_files > 0:
        sys.exit("Aborting due to missing files")

#######################################################################################################################

    ############
    # RUN CODE #
    ############

    pep1 = read_text(args.pep1_filename).rstrip() + "\n"    # Ensures file ends with newline
    pep2 = read_text(args.pep2_filename).rstrip() + "\n"    # Ensures file ends with newline

    if not args.permutations == None:
        print("Building permutation generator for %s permutations..." % (args.permutations))
        permutations = run_permutations(args.permutations, pep1, pep2, args.length, args.specificity)
        print("Built permutation generator!")
        for i in permutations:
            print(len(i.combos))
    elif args.permutations == None:
        concat_alignment = pep1 + pep2
        aln = run_clustal(concat_alignment.encode())
        print(aln)

    sys.exit("Exited successfully!")

#######################################################################################################################