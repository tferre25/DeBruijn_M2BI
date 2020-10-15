#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import matplotlib.pyplot as plt

__author__ = "Théo Ferreira"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Théo Ferreira"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Théo Ferreira"
__email__ = "theo.ferreira.med@gmail.com"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()

#==============================================================
# Création du graphe de Bruijn
#==============================================================

#           Identification des k-mer uniques

def read_fastq(fastq_file) :
    with open (fastq_file) as fillin : 
        for line in fillin :
            yield next(fillin).strip()
            next(fillin)
            next(fillin)


def cut_kmer(sequence, kmer_size) :
    for nt in range(len(sequence) - kmer_size + 1) : 
        yield sequence[nt:nt+kmer_size]


def build_kmer_dict(fastq_file, kmer_size) : 
    kmer_dict = {}
    for sequence in read_fastq(fastq_file):
        for kmer in cut_kmer(sequence, kmer_size):
            if not kmer in kmer_dict : 
                kmer_dict[kmer] = 1
            else : 
                kmer_dict[kmer] += 1      
    return kmer_dict


#           Construction de l'arbre de Bruijn

def build_graph(kmer_dict) :
    G = nx.DiGraph()
    for kmer, weight in kmer_dict.items() : 
        preffix = kmer[:-1]
        suffix = kmer[1:]
        G.add_edge(preffix, suffix ,weight=weight)
    #nx.draw(G, with_labels = True)
    #plt.show()
    return G


#           Parcours du graph 

def get_starting_nodes(G) :
    starting_nodes_list = []
    for node in G.nodes :
        if len(list(G.predecessors(node))) ==  0 :
            starting_nodes_list.append(node)
    return starting_nodes_list

  
def get_sink_nodes(G) :
    sink_nodes_list = []
    for node in G.nodes :
        if len(list(G.successors(node))) ==  0 :
            sink_nodes_list.append(node)
    return sink_nodes_list


def get_contigs(G, starting_nodes_list, sink_nodes_list) : 
    contig_list = []
    for starting_node in starting_nodes_list :
        for sink_node in sink_nodes_list : 
            for path in nx.all_simple_paths(G, 
                source = starting_node, 
                target = sink_node) : 
                contig_list.append((path, len(path)))
    return contig_list

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(contig_list, file_out : str) :
    with open(file_out,'w') as fillout : 
        for i in range(len(contig_list)) : 
            contig = "".join(contig_list[i][0][0])
            for read in range(1, len(contig_list[i][0])) : 
                contig = contig + contig_list[i][0][read][-1]

            fillout.write('>Contig number {} len = {}\n{}\n'
                .format(i, len(contig), fill(contig)))


#==============================================================
# Simplification du graphe de Bruijn
#==============================================================       

#           Résolution des bulles 

def std(values_list): 
    return statistics.stdev(values_list)


def path_average_weight(G, path) :
    weight = 0
    for i in range(len(path) - 1) : 
        weight += graph.edges(path[i], path[i+1], ["weight"])
    return (weight / (len(path) - 1))


def remove_paths(G, path_list) : 
    pass
def delete_entry_node() :
    pass
def select_best_path() : 
    pass
def solve_bubble() :
    pass
def simplify_bubbles() : 
    pass
def solve_entry_tips() : 
    pass
def solve_out_tips() : 
    pass

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

if __name__ == '__main__':
    g = build_graph(build_kmer_dict("../data/eva71_two_reads.fq", 9))
    enter = get_starting_nodes(g)
    exit = get_sink_nodes(g)
    contig = get_contigs(g, enter, exit)
    save_contigs(contig, "output.fasta")
    

    main()
