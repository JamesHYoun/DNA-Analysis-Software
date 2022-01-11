import sys
import argparse
import numpy as np
import time
import zipfile
from collections import defaultdict

def parse_reads_file(reads_fn):
    """
    :param reads_fn: the file containing all of the reads
    :return: outputs a list of all paired-end reads
    """
    try:
        with open(reads_fn, 'r') as rFile:
            print("Parsing Reads")
            first_line = True
            count = 0
            all_reads = []
            for line in rFile:
                count += 1
                if count % 1000 == 0:
                    print(count, " reads done")
                if first_line:
                    first_line = False
                    continue
                ends = line.strip().split(',')
                all_reads.append(ends)
                if len(ends[0]) != 50 or len(ends[1]) != 50:
                    print('error')
        return all_reads
    except IOError:
        print("Could not read file: ", reads_fn)
        return None


def parse_ref_file(ref_fn):
    """
    :param ref_fn: the file containing the reference genome
    :return: a string containing the reference genome
    """
    try:
        with open(ref_fn, 'r') as gFile:
            print("Parsing Ref")
            first_line = True
            ref_genome = ''
            for line in gFile:
                if first_line:
                    first_line = False
                    continue
                ref_genome += line.strip()
        return ref_genome
    except IOError:
        print("Could not read file: ", ref_fn)
        return None


MIN_GAP = 90
MAX_GAP = 110
LEN_READ = 50

OP_PENALTY = -7
EX_PENALTY = -5
MS_PENALTY = -1

START = 'START'
ACROSS = 'ACROSS'
INSERTION = 'INSERTION'
DELETION = 'DELETION'


def count_mismatches(reference, pos, read_part):
    count = 0
    for i in range(len(read_part)):
        if reference[pos + i] != read_part[i]:
            count += 1
    return count


def get_valid_pos(reference, start_pos, read_part):
    for i in range(MAX_GAP - MIN_GAP + 1):
        pos = start_pos + i
        if count_mismatches(reference, pos, read_part) <= 1:
            return pos
    return None


def update_coverage(read, pos_0, pos_1, coverage, success):
    for i in range(len(read[0])):
        coverage[pos_0 + i][read[0][i]] += 1
        coverage[pos_1 + i][read[1][i]] += 1
        if success:
            coverage[pos_0 + i][read[0][i]] = float('inf')
            coverage[pos_1 + i][read[1][i]] = float('inf')
    

def get_full_error_0(indices_a, indices_b):
    for index_a in indices_a:
        for index_b in indices_b:
            diff = index_b - index_a
            if LEN_READ + MIN_GAP <= diff and diff <= LEN_READ + MAX_GAP:
                return True, index_a, index_b
    return False, -1, -1


def get_half_error_0(indices_a, indices_b):
    indices = []
    for index_a in indices_a:
        for index_b in indices_b:
            if index_b - index_a == 25:
                indices.append(index_a)
    return indices


def get_full_error_1(indices_a_e0, indices_b_e1, indices_a_e1, indices_b_e0):
    for index_a_e0 in indices_a_e0:
        for index_b_e1 in indices_b_e1:
            diff = index_b_e1 - index_a_e0
            if LEN_READ + MIN_GAP <= diff and diff <= LEN_READ + MAX_GAP:
                return True, index_a_e0, index_b_e1
    for index_a_e1 in indices_a_e1:
        for index_b_e0 in indices_b_e0:
            diff = index_b_e0 - index_a_e1
            if LEN_READ + MIN_GAP <= diff and diff <= LEN_READ + MAX_GAP:
                return True, index_a_e1, index_b_e0
    return False, -1, -1


def get_half_error_1(reference, read_a, read_b, indices_a, indices_b):
    indices = []
    for index_a in indices_a:
        if count_mismatches(reference, index_a + LEN_READ / 2, read_b) == 1:
            indices.append(index_a)
    for index_b in indices_b:
        if count_mismatches(reference, index_b - LEN_READ / 2, read_a) == 1:
            indices.append(index_b - LEN_READ / 2)
    return indices


def get_full_error_2(indices_a, indices_b):
    for index_a in indices_a:
        for index_b in indices_b:
            diff = index_b - index_a
            if LEN_READ + MIN_GAP <= diff and diff <= LEN_READ + MAX_GAP:
                return True, index_a, index_b
    return False, -1, -1


# Longest Common Subsequence
def lcs_backtrack(v, w):
    n, m = len(v), len(w)
    lower = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    middle = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    upper = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    backtrack = [[START for _ in range(m + 1)] for _ in range(n + 1)]
    for i in range(1, n + 1):
        upper[i][0] = float('-inf')
    for j in range(1, m + 1):
        penalty = EX_PENALTY
        if j == 1:
            penalty = OP_PENALTY
        lower[0][j] = float('-inf')
        middle[0][j] = middle[0][j - 1] + penalty
        backtrack[0][j] = DELETION
    best_score, best_pos = float('-inf'), (-1, -1)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # update lower
            if lower[i-1][j] + EX_PENALTY > middle[i-1][j] + OP_PENALTY:
                lower[i][j] = lower[i-1][j] + EX_PENALTY
            else:
                lower[i][j] = middle[i-1][j] + OP_PENALTY
            # update upper
            if upper[i][j-1] + EX_PENALTY > middle[i][j-1] + OP_PENALTY:
                upper[i][j] = upper[i][j-1] + EX_PENALTY
            else:
                upper[i][j] = middle[i][j-1] + OP_PENALTY
            # update middle
            if lower[i][j] > upper[i][j]:
                middle[i][j] = lower[i][j]
                backtrack[i][j] = INSERTION
            else:
                middle[i][j] = upper[i][j]
                backtrack[i][j] = DELETION
            penalty = MS_PENALTY
            if v[i - 1] == w[j - 1]:
                penalty = 1
            if middle[i-1][j-1] + penalty > middle[i][j]:
                middle[i][j] = middle[i-1][j-1] + penalty
                backtrack[i][j] = ACROSS
            if backtrack[i][j] == START:
                backtrack[i][j] = ACROSS
            if j == m:
                if middle[i][j] > best_score:
                    best_score = middle[i][j]
                    best_pos = (i, j)         
    return middle[best_pos[0]][best_pos[1]], best_pos, backtrack


def get_alignment(v, w, end_i, end_j, backtrack):
    v_n, w_n = '', ''
    i, j = end_i, end_j
    while j > 0:
        if backtrack[i][j] == INSERTION:
            v_n += v[i - 1]
            w_n += '-'
            i -= 1
        elif backtrack[i][j] == DELETION:
            v_n += '-'
            w_n += w[j - 1]
            j -= 1
        elif backtrack[i][j] == ACROSS:
            v_n += v[i - 1]
            w_n += w[j - 1]
            i -= 1
            j -= 1
    return i, (v_n[::-1], w_n[::-1])


def get_indel(reference, indices_a, indices_b, read_a, read_b):
    best_score, best_pos, best_alignment = float('-inf'), -1, ('', '')
    for idx in indices_a:
        start = idx+LEN_READ+MIN_GAP
        score, pos, backtrack = lcs_backtrack(reference[start:start+100], read_b)
        i, j = pos[0], pos[1]
        pos, alignment = get_alignment(reference[start:start+100], read_b, i, j, backtrack)
        if score > best_score and pos <= 20:
            best_score = score
            best_pos = start + pos
            best_alignment = alignment
    for idx in indices_b:
        end = idx-MIN_GAP
        score, pos, backtrack = lcs_backtrack(reference[end-100:end][::-1], read_a[::-1])
        i, j = pos[0], pos[1]
        pos, alignment = get_alignment(reference[end-100:end][::-1], read_a[::-1], i, j, backtrack)
        if score > best_score and pos <= 20:
            best_score = score
            best_pos = end - pos - len(alignment[0])
            best_alignment = (alignment[0][::-1], alignment[1][::-1])
    return best_score, best_pos, best_alignment         
                

def get_coverage(reference, input_reads, mapping):
    count = 0
    coverage = [defaultdict(int) for _ in range(len(reference))]
    insertions, deletions = [], []
    for input_read in input_reads:
        # >SNPS
        
        # first rotation error 0
        read_a = [input_read[0], input_read[1][::-1]]
        indices_a00 = mapping[read_a[0][:25]]
        indices_a01 = mapping[read_a[0][25:]]
        indices_a10 = mapping[read_a[1][:25]]
        indices_a11 = mapping[read_a[1][25:]]
        indices_a0_e0 = get_half_error_0(indices_a00, indices_a01)
        indices_a1_e0 = get_half_error_0(indices_a10, indices_a11)
        status_a_e0, pos_a0_e0, pos_a1_e0 = get_full_error_0(indices_a0_e0, indices_a1_e0)
        if status_a_e0 == True:
            update_coverage(read_a, pos_a0_e0, pos_a1_e0, coverage, True)
            continue
        # second rotation error 0
        read_b = [input_read[0][::-1], input_read[1]]
        indices_b00 = mapping[read_b[0][:25]]
        indices_b01 = mapping[read_b[0][25:]]
        indices_b10 = mapping[read_b[1][:25]]
        indices_b11 = mapping[read_b[1][25:]]
        indices_b0_e0 = get_half_error_0(indices_b00, indices_b01)
        indices_b1_e0 = get_half_error_0(indices_b10, indices_b11)
        status_b_e0, pos_b0_e0, pos_b1_e0 = get_full_error_0(indices_b0_e0, indices_b1_e0)
        if status_b_e0 == True:
            update_coverage(read_b, pos_b0_e0, pos_b1_e0, coverage, True)
            continue
        # first rotation error 1
        indices_a0_e1 = get_half_error_1(reference, read_a[0][:25], read_a[0][25:], indices_a00, indices_a01)
        indices_a1_e1 = get_half_error_1(reference, read_a[1][:25], read_a[1][25:], indices_a10, indices_a11)
        status_a_e1, pos_a0_e1, pos_a1_e1 = get_full_error_1(indices_a0_e1, indices_a1_e0, indices_a0_e0, indices_a1_e1)
        if status_a_e1 == True:
            update_coverage(read_a, pos_a0_e1, pos_a1_e1, coverage, False)
            continue
        # second rotation error 1
        indices_b0_e1 = get_half_error_1(reference, read_b[0][:25], read_b[0][25:], indices_b00, indices_b01)
        indices_b1_e1 = get_half_error_1(reference, read_b[1][:25], read_b[1][25:], indices_b10, indices_b11)
        status_b_e1, pos_b0_e1, pos_b1_e1 = get_full_error_1(indices_b0_e1, indices_b1_e0, indices_b0_e0, indices_b1_e1)
        if status_b_e1 == True:
            update_coverage(read_b, pos_b0_e1, pos_b1_e1, coverage, False)
            continue
        # first rotation error 2
        status_a_e2, pos_a0_e2, pos_a1_e2 = get_full_error_2(indices_a0_e1, indices_a1_e1)
        if status_a_e2 == True:
            update_coverage(read_a, pos_a0_e2, pos_a1_e2, coverage, False)
            continue
        # second rotation error 2
        status_b_e2, pos_b0_e2, pos_b1_e2 = get_full_error_2(indices_b0_e1, indices_b1_e1)
        if status_b_e2 == True:
            update_coverage(read_b, pos_b0_e2, pos_b1_e2, coverage, False)
            continue
        
        # >INDELS
        count += 1
        best_score, best_pos, best_alignment = float('-inf'), -1, ('','')
        # first rotation
        score, pos, alignment = get_indel(reference, indices_a0_e0, indices_a1_e0, read_a[0], read_a[1])
        if score > best_score:
            best_score = score
            best_pos = pos
            best_alignment = alignment
        # second rotation
        score, pos, alignment = get_indel(reference, indices_b0_e0, indices_b1_e0, read_b[0], read_b[1])
        if score > best_score:
            best_score = score
            best_pos = pos
            best_alignment = alignment
             
        ins_start, del_start = False, False
        pos_ref = best_pos
        z = ''
        for i in range(0, len(best_alignment[0])):
            if best_alignment[0][i] != '-' and ins_start:
                ins_start = False
                insertions.append([z, pos_ref])
                z = ''
            if best_alignment[1][i] != '-' and del_start:
                del_start = False
                deletions.append([z, pos_ref - len(z)])
                z = ''
            if best_alignment[0][i] == '-':
                ins_start = True
                z += best_alignment[1][i]
            if best_alignment[1][i] == '-':
                del_start = True
                z += best_alignment[0][i]
            if not ins_start:
                pos_ref += 1
        if ins_start:
            insertions.append([z, pos_ref])
        if del_start:
            deletions.append([z, pos_ref - len(z)])      
    return coverage, insertions, deletions


def get_snps(reference, coverage):
    snps = []
    for i in range(len(reference)):
        base_counts = coverage[i]
        max_base, max_count = reference[i], float('-inf')
        for base in base_counts:
            if base_counts[base] > max_count:
                max_base = base
                max_count = base_counts[base]
        if reference[i] != max_base:
            snps.append([reference[i], max_base, i])
    return snps


def get_mapping(input_file):
    mapping = defaultdict(list)
    with open(input_file, 'r') as f:
        line = f.readline()[:-1]
        while line != '':
            line = line.split(' ')
            mapping[line[0]] += [int(x) for x in line[1:]]
            line = f.readline()[:-1]
    return mapping


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_hasher.py takes in data consisting '
                                     'of a genome and a set of reads and aligns the reads to the reference genome, '
                                     'then calls SNPS and indels based on this alignment.')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) practice_W_3_chr_1 for 10K length genome practice data\n'
                             '2) practice_E_1_chr_1 for 1 million length genome practice data\n'
                             '3) hw2undergrad_E_2_chr_1 for project 2 undergrad for-credit data\n'
                             '4) hw2grad_M_1_chr_1 for project 2 grad for-credit data\n')
    args = parser.parse_args()
    reference_fn = args.reference_file
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)
    reference = parse_ref_file(reference_fn)
    if reference is None:
        sys.exit(1)

    mapping_file = ''
    if args.output_header == 'practice_W_3_chr_1':
        mapping_file = 'mapping_W_3.txt'
    if args.output_header == 'practice_E_1_chr_1':
        mapping_file = 'mapping_E_1.txt'
    if args.output_header == 'hw2undergrad_E_2_chr_1':
        mapping_file = 'mapping_E_2.txt'

    mapping = get_mapping(mapping_file)
    coverage, insertions, deletions = get_coverage(reference, input_reads, mapping)
    insertions.sort(key=lambda x:x[1])
    deletions.sort(key=lambda x:x[1])
    snps = get_snps(reference, coverage)

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n>SNP\n')
        for x in snps:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>INS\n')
        for x in insertions:
            output_file.write(','.join([str(u) for u in x]) + '\n')
        output_file.write('>DEL\n')
        for x in deletions:
            output_file.write(','.join([str(u) for u in x]) + '\n')
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
