# -*- coding: utf-8 -*-
"""
logol_analyse provide some analyse tools for logol xml results.
Without any option, it will provide the number of hit, how many sequences have
at least one hit, and a graph with the repartition of the hits.

Usage:
    logol_analyse.py <input> <data> [options]

options:
    --graph, -g=<name>       The graph name, to save it directly.
    --help, -h               It call help. UNBELIEVABLE!!!!!
    --nograph -n             No graph creation
    --origin, -o INT         The 0 emplacement on sequences [default: 150]
    --position -p=<name>     Return a file containing position of each motif
    --result -r=<name>       Save a fasta file with the matched sequences.
    --signature, -s=<name>   Create a file with for each sequences the hits.
    --hits, -t               Display a hits/sequences graph.
    --version, -v            Maybe it's a trap ^^
    --xclude, -x=<name>      Create a file containing all unmatched sequences

"""

##########
# IMPORT #
##########
import matplotlib.pyplot as plt
import pylab
import glob
import os

from docopt import docopt
from lxml import etree
from Bio import SeqIO

#############
# ARGUMENTS #
#############
if __name__ == '__main__':
    arguments = docopt(__doc__, version = '1.3')

########
# MAIN #
########
def __main__(arguments):
    total        = 0
    count        = 0
    hit          = []
    se           = set()  # Contain sequences header
    hits_per_seq = []
    # Here we check all the .xml file
    for f in glob.glob(os.getcwd()+"/"+arguments['<input>']+"*.xml"):
        nb_hit = 0
        total += 1
        tree = etree.parse(f)
        # Collect of the hit beginning and ID
        for seq in tree.xpath("/sequences/match/begin"):
            count += 1
            nb_hit +=1
            hit.append(int(seq.text)-int(arguments['--origin']))
            [se.add(a.text) for a in tree.xpath("/sequences/fastaHeader")]
        if nb_hit > 0:
            hits_per_seq.append(nb_hit)
    print("Nombre de hits: "+str(count))
    print("Nombre de séquences touchées: "+str(len(se))+" sur "+str(total))
    print("Nombre max de hits par séquences: "+str(max(hits_per_seq)))
    if arguments['--result'] != None:
        seq_match(se)
    if arguments['--xclude'] != None:
        seq_no_match(se)
    if arguments['--nograph'] == False:
        graph(hit)
    if arguments['--signature'] != None:
        save_signature()
    if arguments['--position'] != None:
        save_position()
    if arguments['--hits'] != False:
        display_hits(hits_per_seq)

#############
# FUNCTIONS #
#############
def seq_match(seq):
    out  = open(os.getcwd()+'/'+arguments['--result'], 'w')
    data = open(os.getcwd()+'/'+arguments['<data>'], "rU")
    for s in SeqIO.parse(data, "fasta"):
        if s.id in seq:
            out.write(s.format("fasta"))
    out.close()
    data.close()

def seq_no_match(seq):
    out  = open(os.getcwd()+'/'+arguments['--xclude'], 'w')
    data = open(os.getcwd()+'/'+arguments['<data>'], "rU")
    for s in SeqIO.parse(data, "fasta"):
        if s.id not in seq:
            out.write(s.format("fasta"))
    out.close()
    data.close()

def graph(hit):
    plt.hist(hit, range(min(hit), max(hit)))
    plt.xticks(range(min(hit), max(hit), 10))
    plt.xlabel("Emplacement des hits sur les séquences")
    plt.ylabel("Nombre de hits")
    if arguments['--graph'] != None:
        plt.savefig(arguments['--graph']+'.png')
        pylab.close()
    else:
        plt.show()

def save_signature():
    sign = open(os.getcwd()+'/'+arguments['--signature'], 'w')
    for f in glob.glob(os.getcwd()+"/"+arguments['<input>']+"*"):
        fr   = []  # Will have the last char of var, which is frag nb
        c    = 0
        tree = etree.parse(f)
        if  tree.xpath("/sequences/match/variable") != []:
            [sign.write('>'+h.text+'\n') for h in tree.xpath("/sequences/fastaHeader")]
            [fr.append((int(i.get("name")[-1]))) for i in tree.xpath("/sequences/match/variable")]
            m = max(fr) # Fragments number to have the complete match
            for i in tree.xpath("/sequences/match/variable/content"):
                c += 1
                sign.write(i.text)
                if c >= m:
                    sign.write("\n")
                    c = 0
    sign.close()

def save_position():
    begin = [] # Will contain all the begining number
    end   = []
    seq   = [] # Will contain all the sequences found
    iD    = [] # Will contair the sequences ID
    n     = 0  # nb of line we will have to write
    i     = 0
    pos   = open(os.getcwd()+'/'+arguments['--position'], 'w')
    pos.write("ID\tbegin\tsequence\tend\n") 
    for f in glob.glob(os.getcwd()+"/"+arguments['<input>']+"*"):
        tree = etree.parse(f)
        for s in tree.xpath("/sequences/match/variable/content"):
            n += 1
            seq.append(s.text)
            [iD.append(h.text) for h in tree.xpath("/sequences/fastaHeader")]
        for b in tree.xpath("/sequences/match/variable/begin"):
            begin.append(str(b.text))
        for e in tree.xpath("/sequences/match/variable/end"):
            end.append(str(e.text))
    # Now, we write those info into the file
    while i < n:
        pos.write(iD[i]+"\t"+begin[i]+"\t"+seq[i]+"\t"+end[i]+"\n")
        i += 1
    pos.close()

def display_hits(hits_per_seq):
    plt.hist(hits_per_seq, range(min(hits_per_seq), max(hits_per_seq)))
    plt.xticks(range(min(hits_per_seq), max(hits_per_seq), 1))
    plt.xlabel("Nombre de hits par séquences")
    plt.ylabel("Nombre de séquences")
    plt.show()

##########
# LAUNCH #
##########
__main__(arguments)
