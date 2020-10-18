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


import statistics
import argparse
import os
import sys
import random
import networkx as nx
random.seed(9001)

__author__ = "Théo Ferreira"
__copyright__ = "Universite de Paris"
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
    parser = argparse.ArgumentParser(description=__doc__, usage="{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()

# ==============================================================
# Création du graph de Bruijn
# ==============================================================

#           Identification des k-mer uniques

def read_fastq(fastq_file):
    '''
    Function to read fastq files
    args :
        fastq_file : file with fastq extension
    returns a generator of the sequences
    '''
    with open(fastq_file) as fillin:
        for line in fillin:
            yield next(fillin).strip()
            next(fillin)
            next(fillin)

def cut_kmer(sequence, kmer_size):
    '''
    Function to cut the sequence in kmer with a k size
    args :
        sequence : sequences in fastq
        kmer_size : size of kmer, default 21
    returns a generator of kmers of the sequence
    '''
    for letter in range(len(sequence) - kmer_size + 1):
        yield sequence[letter: letter + kmer_size]

def build_kmer_dict(fastq_file, kmer_size):
    '''
    Function to allocate a kmer to its occurence in a dictionnary
    args :
        fastq_file : file with fastq extension
        kmer_size : size of kmer, default 21
    returns a dictionnary {key = kmer : values = occurence}
    '''
    kmer_dict = {}
    for sequence in read_fastq(fastq_file):
        for kmer in cut_kmer(sequence, kmer_size):
            if not kmer in kmer_dict:
                kmer_dict[kmer] = 1
            else:
                kmer_dict[kmer] += 1
    return kmer_dict

#           Construction de l'arbre de Bruijn

def build_graph(kmer_dict):
    '''
    Fucntion to build a directed graph of suffix and preffix kmer
    args :
        kmer_dict : a dictionnary with kmers and their occurence
    returns a directed graph from kmer preffix to kmer suffix
    '''
    graph = nx.DiGraph()
    for kmer, weight in kmer_dict.items():
        preffix = kmer[:-1]
        suffix = kmer[1:]
        graph.add_edge(preffix, suffix, weight=weight)
    return graph

#           Parcours du graph

def get_starting_nodes(graph):
    '''
    Function to get entry nodes
    args :
        graph : graph
    returns a list of entry nodes
    '''
    starting_nodes_list = []
    for node in graph.nodes:
        if len(list(graph.predecessors(node))) == 0:
            starting_nodes_list.append(node)
    return starting_nodes_list

def get_sink_nodes(graph):
    '''
    Function to get exit nodes
    args :
        graph : graph
    returns a list of exit nodes
    '''
    sink_nodes_list = []
    for node in graph.nodes:
        if len(list(graph.successors(node))) == 0:
            sink_nodes_list.append(node)
    return sink_nodes_list

def get_contigs(graph, starting_nodes_list, sink_nodes_list):
    '''
    Fucntion to give possible paths between entry nodes and exit nodes
    args :
        G : graph
        starting_nodes_list : list of entry nodes
        sink_nodes_list : list of exit nodes
    returns a list of tuples (contig, length of contig)
    '''
    contig_list = []
    for starting_node in starting_nodes_list:
        for sink_node in sink_nodes_list:
            for path in nx.all_simple_paths(graph,
                                            source=starting_node,
                                            target=sink_node):
                contig = path[0]
                for i in range(1, len(path)):
                    contig += path[i][-1]

                contig_list.append((contig, len(contig)))
    return contig_list

def fill(text, width=80):
    '''
    Split text with a line return to respect fasta format
    '''
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(contig_list, file_out: str):
    '''
    Function to save contigs in txt file (fasta format)
    args :
        contig_list : list of tuples (contig, length of contig)
        file_out : name of the output file
    returns a file (fasta format)
    '''
    with open(file_out, 'w') as fillout:
        for i in range(len(contig_list)):
            contig = contig_list[i][0]

            fillout.write('>contig_{} len={}\n{}\n'.format(
                i, len(contig), fill(contig)))

# ==============================================================
# Simplification du graphe de Bruijn
# ==============================================================

#           Résolution des bulles

def std(values_list):
    '''
    returns the standard deviation
    '''
    return statistics.stdev(values_list)

def path_average_weight(graph, path):
    '''
    Function to get the mean weight of a path (average occurrence)
    args :
        graph : graph
        path : a path in the graph (possible contig)
    returns the average weight of edges in the path (arg)
    '''
    weight = 0
    for node in range(len(path) - 1):
        weight += graph.edges[path[node], path[node + 1]]["weight"]
    return weight / (len(path) - 1)

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    '''
    Funtion to remove all path from a path list
    args :
        graph : graph
        path_list : Path list
        delete_entry_node: Entry node deletion boolean
        delete_sink_node: Sink node deletion boolean
    returns graph
    '''
    for node in range(len(path_list)):
        # supprime noeuds entre entrée et sortie
        graph.remove_nodes_from(path_list[node][1:-1])
        if delete_entry_node == True:
            graph.remove_node(path_list[node][0])
        if delete_sink_node == True:
            graph.remove_node(path_list[node][-1])
    return graph

def select_best_path(graph, path_list,
                     path_len_list, path_weight_list,
                     delete_entry_node=False, delete_sink_node=False):
    '''
    Function to select the best path (maximun weight, then maximum length)
    if there was more than one better path, the path was chosen at random
    args :
        graph : graph
        path_list : Path list
        path_len_list : path length list
        path_weight_list : average path weight list
        delete_entry_node: Entry node deletion boolean
        delete_sink_node: Sink node deletion boolean
    returns a grah with the best path
    '''

    # Weight
    weight_idx = []
    best_paths = []

    for path, weight in enumerate(path_weight_list):
        if weight == max(path_weight_list):
            weight_idx.append(path)

    # Adding a condition : if len > 1
    # because variance can only be calculated on less than two values
    if len(path_weight_list) > 1:
        if std(path_weight_list) > 0:
            best_paths = [path_list[i] for i in weight_idx]
        else:
            best_paths = path_list
    else:
        best_paths = path_list

    # Length
    if len(path_len_list) > 1:
        if std(path_len_list) > 0:
            path_len_list = [path_len_list[i] for i in weight_idx]
            len_idx = []
            for path, length in enumerate(path_len_list):
                if length == max(path_len_list):
                    len_idx.append(path)
            best_paths = [best_paths[i] for i in len_idx]
        else:
            best_paths = best_paths
    else:
        best_paths = best_paths

    # Random    
    if len(best_paths) > 1:
        best_paths = best_paths[random.randint(0, len(best_paths))]
    path_list.remove(best_paths[0])

    graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    return graph

def solve_bubble(graph, node_ancestor, node_suc):
    '''
    Function to select the best path in a bubble
    args :
        garph : graph
        node_ancestor : Ancestor node
        node_suc : Successor node
    returns graph with best path in a bubble
    '''

    bubble_path_list   = \
    list(nx.all_simple_paths(graph,node_ancestor,node_suc))
    bubble_path_len    = \
    [len(path) for path in nx.all_simple_paths(graph,node_ancestor,node_suc)]
    bubble_path_weight = \
    [path_average_weight(graph, path) for path in nx.all_simple_paths(graph,node_ancestor,node_suc)]

    graph = select_best_path(
        graph,
        bubble_path_list,
        bubble_path_len, bubble_path_weight
    )
    return graph

def simplify_bubbles(graph):
    '''
    Funtion to remove all bubbles
    args :
        graph : graph
    returns graph with no bubble
    '''
    graph_simplify = graph.copy()
    for node in graph.nodes:
        predecessors_list = list(graph.predecessors(node))
        while len(predecessors_list) > 1:
            common_ancestor = nx.lowest_common_ancestor(graph,
                                                        predecessors_list[0],
                                                        predecessors_list[1]
                                                        )
            graph_simplify = solve_bubble(graph_simplify, common_ancestor, node)
            predecessors_list = list(graph_simplify.predecessors(node))
    return graph_simplify

def solve_entry_tips(graph, starting_nodes_list):
    '''
    Funtion to remove entry tips
    args :
        graph : graph
        starting_nodes_list : list of entry nodes
    returns a graph with no entry tips
    '''
    path_weight_list = []
    path_len_list = []
    path_list = []

    for node_source in starting_nodes_list:
        if len(list(graph.successors(node_source))) > 0:
            node_succ = list(graph.successors(node_source))[0]
            nodes_pred = list(graph.predecessors(node_succ))

            while len(nodes_pred) == 1 \
                    and len(list(graph.successors(node_succ))) > 0:
                node_succ  = list(graph.successors(node_succ))[0]
                nodes_pred = list(graph.predecessors(node_succ))

            for path in nx.all_simple_paths(graph, node_source, node_succ):
                path_list.append(path)
                path_weight_list.append(path_average_weight(graph, path))
                path_len_list.append(len(path))

    graph = select_best_path(
        graph,
        path_list,
        path_weight_list,
        path_len_list,
        delete_entry_node=True
    )
    return graph

def solve_out_tips(graph, sink_nodes_list):
    '''
    Funtion to remove out tips
    args :
        graph : graph
        starting_nodes_list : list of out nodes
    returns a graph with no out tips
    '''
    path_weight_list = []
    path_len_list = []
    path_list = []

    for node_target in sink_nodes_list:
        if len(list(graph.predecessors(node_target))) > 0:
            node_pred = list(graph.predecessors(node_target))[0]
            nodes_succ = list(graph.successors(node_pred))

            while len(nodes_succ) == 1 \
                    and len(list(graph.predecessors(node_pred))) > 0:
                node_pred = list(graph.predecessors(node_pred))[0]
                nodes_succ = list(graph.successors(node_pred))

            for path in nx.all_simple_paths(graph, node_pred, node_target):
                path_list.append(path)
                path_weight_list.append(path_average_weight(graph, path))
                path_len_list.append(len(path))

    graph = select_best_path(
        graph,
        path_list,
        path_weight_list,
        path_len_list,
        delete_sink_node=True
    )
    return graph


# ==============================================================
# Main program
# ==============================================================

def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Get file
    read_fastq(args.fastq_file)

    # Build kmer dictionnary
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)

    # Build directed graph
    graph = build_graph(kmer_dict)

    # Simplify the graph
    graph = simplify_bubbles(graph)

    # Solve bubbles : One entry - One sink, remove non optimal nodes
    starting_nodes = get_starting_nodes(graph)
    sink_nodes = get_sink_nodes(graph)
    graph = solve_entry_tips(graph, starting_nodes)
    graph = solve_out_tips(graph, sink_nodes)

    # Get contigs and save it
    starting_nodes = get_starting_nodes(graph)
    sink_nodes = get_sink_nodes(graph)
    contigs = get_contigs(graph, starting_nodes, sink_nodes)
    save_contigs(contigs, args.output_file)


if __name__ == '__main__':
    main()
