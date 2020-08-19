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
    if len(homology_regions) == 0:
        return [(0, alignment_length)]
    else:
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

def run_alignment(n_permutations, pep1, pep2, length, specificity):
    if n_permutations == 0:
        yield clustal_alignment(0, pep1, pep2, length, specificity)
    elif n_permutations > 0:
        for _ in range(n_permutations):
            yield clustal_alignment(random.randint(1,2), pep1, pep2, length, specificity)

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

class clustal_alignment:
    help = 'clustal alignment of permuted protein sequence'
    # zero for unpermuted, 1 for permute seq 1, 2 for permute seq 2
    def __init__(self, permute, pep1, pep2, length, specificity): # 0 for unpermuted, 1 for permute seq 1
        if permute == 0:
            self.pep1 = pep1
            self.pep2 = pep2
        elif permute == 1:
            self.pep1 = pep1
            self.pep2 = permute_peptide_fasta(pep2)
        elif permute == 2:
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
        self.non_homology_combos = list(return_combos(self.offset1, self.offset2, self.non_homology_regions))
        self.homology_combos = list(return_combos(self.offset1, self.offset2, self.homology_regions))

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
    parser.add_argument("-A", "--all", 
                        help='''By default, CHIMERA will generate a list for only
                                potential chimeras where transitions between homologs
                                occurs at a point of confident homology. Setting the
                                --all flag will also generate chimeras with transitions
                                outside of regions of confident homology. Used with
                                -n, --indel. See <figure> for explanation.''',
                        action="store_true")
    parser.add_argument("-n", "--indel",
                        type=int,
                        nargs='?',
                        const=1,
                        default=0,
                        help='''Defines maximum number of added/removed
                                codons, up to a (but not including) the next region
                                of confident homology. See <figure> for explanation.
                                To include all possible transitions, set to -1.
                                 Default: 0''')
    parser.add_argument("-L", "--length",
                        type=int,
                        nargs='?',
                        const=1,
                        default=4,
                        help='''Minimum number of consecutive amino acids above the
                                specificity threshold to determine a region of
                                confident homology (defined by -S). Default: 4''')
    parser.add_argument("-s", "--specificity",
                        type=str,
                        nargs='?',
                        const=1,
                        default="high",
                        help='''Accepts 'low' 'medium' or 'high'.
                                Defines degree of amino acid similarity from clustal alignment
                                for defining regions of confident homology. Used with -L.
                                'low' includes [.:*], 'medium includes' [:*], 'high' includes [*].
                                Default: high.''')
    parser.add_argument("-p", "--permutations",
                        type=int,
                        nargs='?',
                        const=1,
                        help = '''Number of alignment permutations to run. Default: 100.''')
    parser.add_argument("-u", "--unique", 
                        help='''Accepts 'prot' or 'dna'. 
                                Default 'prot' will deduplicate oligo sequences, retaining one representative
                                for a given protein sequence. 'dna' will retain identical dna sequences,
                                with the possibility of duplicate resulting protein sequences (which may be of
                                interest for codon bias).''',
                        type=str,
                        nargs='?',
                        const=1,
                        default="prot")

    args = parser.parse_args()
    if args.debug:
        debug = True
        print("debug mode on, increasing verbosity.")
    else:
        debug = False

    if debug == True:
        print("""running with length = %s and specificity = %s""" % (args.length, args.specificity))

#######################################################################################################################

    ######################
    # VALIDATE ARGUMENTS #
    ######################

    try: assert args.length > 0, "ERROR: length must be > 0"
    except AssertionError as error: sys.exit(error)

    try: assert args.specificity in ["low", "medium", "high"], "ERROR: -s, --specificity must be 'low', 'medium', or 'high'"
    except AssertionError as error: sys.exit(error)

    try: assert args.unique in ["prot", "dna"], "ERROR: -u, --unique must be 'prot' or 'dna'"
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
        permutations = run_alignment(args.permutations, pep1, pep2, args.length, args.specificity)
        print("Built permutation generator!")
        for i in permutations:
            print(i.homology_regions)
            print(i.homology_combos)
        # with open(args.specificity + '.' + str(args.length) + '.permutations', 'w') as outfile:
        #     outfile.write(                           
        #          '\t'.join(["specificity", "match_length", "homology_regions", "alignment_length", "n_chimeras"]) + '\n'
        #                 )
        #     for i in permutations:
        #         outfile.write(
        #                     '\t'.join([str(x) for x in [args.specificity, args.length, len(i.homology_regions), len(i.aln1), len(i.combos)]]) + '\n'
        #                     )
    elif args.permutations == None:
        alignment = run_alignment(0, pep1, pep2, args.length, args.specificity)
        i = list(alignment)[0]
        if not os.path.isfile("true_alignment.out"):
            with open("true_alignment.out", 'w') as outfile:
                outfile.write(                           
                 '\t'.join(["specificity", "match_length", "homology_regions", "alignment_length", "n_chimeras"]) + '\n'
                        )
        with open("true_alignment.out", 'a') as outfile:
            outfile.write('\t'.join([str(x) for x in [args.specificity, args.length, len(i.homology_regions), len(i.aln1), len(i.combos)]]) + '\n')

    sys.exit("Exited successfully!")

#######################################################################################################################
# NOTES / TODO

# Round 1 permutations to determine area to focus in on
# Round 2 permutations to get best FDR
# Option to just get indel = 0 chimeras
# Option to only make chimeras within high confidence homology regions


# Recently done:
# added chimeras within the homology regions class.homology_combos
# add speedup to ignore non_homology_combos if options don't want them