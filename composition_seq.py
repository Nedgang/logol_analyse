# -*- coding: utf-8 -*-
"""
composition_seq analyse a file containing DNA or RNA and return a graph with 
the sequences base composition.

Usage:
    composition_seq <file> [options]

Options:
    --format, -f=<name>  The file format type [default: fasta].
    --help, -h           Display this.
    --rna, -r            If it is RNA sequences.
    --version, -v        Script version.
    --window, -w=<int>   The size of the window [default: 1].

"""

##########
# IMPORT #
##########
import matplotlib.pyplot as plt
import pylab
import os

from docopt import docopt
from Bio import SeqIO

#############
# ARGUMENTS #
#############
if __name__ == '__main__':
    arguments = docopt(__doc__, version = '2.0')

########
# MAIN #
########
def __main__(arguments):
    list_sequences = _FASTA_reader(arguments["<file>"], arguments["--format"])
    scan           = _sequence_scan(list_sequences, int(arguments["--window"]))
    if int(arguments["--window"]) > 1:
        scan       = _adaptation_for_graph(scan, int(arguments["--window"]))
    _draw_graph(scan, os.path.basename(arguments["<file>"]), arguments["--rna"])


#############
# FUNCTIONS #
#############
def _FASTA_reader(path_file, file_format):
    """
    Return a list  of all the sequences contained by the file.
    path_file: the path to the file
    file_format: his format
    """
    print("Extracting the sequences.")
    seq = []
    for sequence in SeqIO.parse(os.path.abspath(path_file), file_format):
        seq.append(sequence.seq)
    return seq

def _sequence_scan(list_sequences, window):
    """
    Return a list of 4 lists,(A, T, C, G), containing for each one the 
    frequence of a nucleotide in each position of the sequences.
    list_sequences: a list of DNA sequences.
    window: the size of the analyse window
    """
    print("Sequences analyse.")
    base_value       = float((1/len(list_sequences)/window))
    qtA              = []
    qtT              = []
    qtC              = []
    qtG              = []
    scan             = [qtA, qtT, qtG, qtC]
    scan_emplacement = 0
    base             = 0  # The nt where I am in the sequence
    # Initialisation of the scan window
    while base < len(max(list_sequences)):
        qtA.append(0)
        qtT.append(0)
        qtG.append(0)
        qtC.append(0)
        # For each sequence long enough
        sequences = [seq for seq in list_sequences if len(seq) > base+window]
        for seq in sequences:
            # Letter analyse in the current window
            for letter in seq[base:base+window]:
                if letter == "A":
                    qtA[scan_emplacement] += base_value
                elif letter == "C":
                    qtC[scan_emplacement] += base_value
                elif letter == "G":
                    qtG[scan_emplacement] += base_value
                else:
                    qtT[scan_emplacement] += base_value
        # mark incrementation
        base             += window
        scan_emplacement += 1
    return scan

def _adaptation_for_graph(scan, window):
    """
    In case of a window > 1, the data provide by _sequence_scan is not
    effectivelly used by _draw_graph and need an adaptation.
    scan: a list of 4 list (A, T, G, C)
    window: size of the window
    """
    print("Data adaptation for graph.")
    adapt    = [[],[],[],[]]
    position = 0
    # For each list in scan
    for l in scan:
        for stat in l:
            i = 0
            while i < window:
                adapt[position].append(stat)
                i += 1
        position += 1
    return adapt

def _draw_graph(scan, file_name, rna):
    """
    Draw a graph with the nucleotide composition of sequences.
    scan: the result of nucleotide analyse, a list of four list, A, T, G, C
    file_name: the name of the sequence file
    rna: if it's RNA sequences
    """
    print("Drawing the graph.")
    # Draw the line for each nt
    plt.plot(scan[0], label="A")
    if rna:
        plt.plot(scan[1], label="U")
    else:
        plt.plot(scan[1], label="T")
    plt.plot(scan[2], label="G")
    plt.plot(scan[3], label="C")
    # Display labels
    plt.legend()
    plt.xlabel("Sequence (nt)")
    plt.ylabel("Proportion")
    plt.title("Nucléotide proportion from "+file_name)
    # The size of x axis
    plt.xticks(range(0, len(scan[0]), 10))
    # Show the graph
    plt.show()

##########
# LAUNCH #
##########
__main__(arguments)
