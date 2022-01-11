from os.path import join
import sys
import time
from collections import defaultdict, Counter
import sys
import os
import zipfile
import argparse
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))


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
        return all_reads
    except IOError:
        print("Could not read file: ", reads_fn)
        return None


# Construct de Bruijn Graph.
def debruijn(patterns):
    bruijn = defaultdict(dict)
    for pattern in patterns:
        if pattern[1:] not in bruijn[pattern[:-1]]:
            bruijn[pattern[:-1]][pattern[1:]] = 1
        else:
            bruijn[pattern[:-1]][pattern[1:]] += 1
    return bruijn


# Construct and utilize Eulerian Cycle.
def dfs(egraph, cycle, visited, curr, start):
    while True:
        for node in egraph[curr]:
            if (curr, node) not in visited:
                cycle.append(curr)
                egraph[curr][node] -= 1
                if egraph[curr][node] == 0:
                    visited.add((curr, node))
                outdeg[curr] -= 1
                if start == node or outdeg[node] == 0:
                    return
                curr = node
                break


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_assembly.py takes in data for homework assignment 3 consisting '
                                                 'of a set of reads and aligns the reads to the reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')
    parser.add_argument('-t', '--outputHeader', required=True, dest='output_header',
                        help='String that needs to be outputted on the first line of the output file so that the\n'
                             'online submission system recognizes which leaderboard this file should be submitted to.\n'
                             'This HAS to be one of the following:\n'
                             '1) spectrum_A_1_chr_1 for 10K spectrum reads practice data\n'
                             '2) practice_A_2_chr_1 for 10k normal reads practice data\n'
                             '3) hw3all_A_3_chr_1 for project 3 for-credit data\n')
    args = parser.parse_args()
    reads_fn = args.reads_file

    input_reads = parse_reads_file(reads_fn)
    if input_reads is None:
        sys.exit(1)

    temp = defaultdict(int)
    total = 0
    k = 26
    for read in input_reads:
        for i in range(len(read[0]) - k + 1):
            temp[read[0][i:i+k]] += 1
            total += 1
        for i in range(len(read[1]) - k + 1):
            temp[read[1][::-1][i:i+k]] += 1
            total += 1
    patterns = []
    for x in temp:
        if temp[x] / total > 0.000015:
            patterns.append(x)
    egraph = debruijn(patterns)
    outdeg = defaultdict(int)
    indeg = defaultdict(int)
    odd_nodes = ["", ""]
    for parent in egraph:
        children = egraph[parent]
        for child in children:
            outdeg[parent] += egraph[parent][child]
            indeg[child] += egraph[parent][child]
    starts = []
    for parent in egraph:
        if outdeg[parent] > indeg[parent]:
            starts.append(parent)
    if len(starts) == 0:
        starts.append(egraph[0])
    print(len(starts))
    
    contigs = []
    visited = set()
    increment_0 = 1
    increment_1 = len(starts) // 128
    for i in range(0, len(starts), increment_0):
        start = starts[i]
        ecycle = []
        dfs(egraph, ecycle, visited, start, start)
        text = ''
        for node in ecycle:
            text += node[0]
        last = ecycle[-1]
        for i in range(1, len(last)):
            text += last[i]
        contigs.append(text)

    output_fn = args.output_file
    zip_fn = output_fn + '.zip'
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + args.output_header + '\n')
        output_file.write('>ASSEMBLY\n')
        output_file.write('\n'.join(contigs))
    with zipfile.ZipFile(zip_fn, 'w') as myzip:
        myzip.write(output_fn)
