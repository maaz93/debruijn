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
from operator import itemgetter
import random
import statistics
import matplotlib
import matplotlib.pyplot as plt
random.seed(9001)
from random import randint
import networkx as nx

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
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
                        default=22, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file) as f:
        for line in f :
            yield next(f).rstrip("\n")
            next(f)
            next(f)

def cut_kmer(read, kmer_size):
    for i in range(len(read) - kmer_size + 1):
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    seq=read_fastq(fastq_file)
    kmer_dict = dict()
    for i in seq :
        for j in cut_kmer(i,kmer_size):
            if j in kmer_dict:
                kmer_dict[j] += 1
            else:
                kmer_dict[j] = 1
    return kmer_dict



def build_graph(kmer_dict):
    G = nx.DiGraph()
    for i in kmer_dict:
        G.add_edge(i[:-1],i[1:],weight=kmer_dict[i])
    return G


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    if delete_entry_node is True and delete_sink_node is True:
        for path in path_list:
            graph.remove_nodes_from(path)
    elif delete_entry_node is True and delete_sink_node is False:
        for path in path_list:
            graph.remove_nodes_from(path[:-1])
    elif delete_entry_node is False and delete_sink_node is True:
        for path in path_list:
            graph.remove_nodes_from(path[1:])
    else:
        for path in path_list:
            graph.remove_nodes_from(path[1:-1])
    return graph



def std(data):
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    path1 = [path for path in path_list]
    if std(weight_avg_list) != 0:
        max_weight = max(weight_avg_list)
        path1.pop(weight_avg_list.index(max_weight))
        return remove_paths(graph,path1,delete_entry_node,delete_sink_node)

    elif std(path_length) != 0:
        max_length = max(path_length)
        path1.pop(path_length.index(max_length))
        return remove_paths(graph,path1,delete_entry_node,delete_sink_node)
    else :
        path1 = path_list.pop(path_list[random.randint()])
        return remove_paths(graph,path1,delete_entry_node,delete_sink_node)

def path_average_weight(graph, path):
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])


def solve_bubble(graph, ancestor_node, descendant_node):
    path_list1 = nx.all_simple_paths(graph,ancestor_node,descendant_node)
    path_list = [path for path in path_list1]
    weight_avg_list = [path_average_weight(graph,path) for path in path_list]
    path_length = [len(path) for path in path_list]
    return select_best_path(graph,path_list,path_length,weight_avg_list)

def simplify_bubbles(graph):
    bubble = False
    for n in graph.nodes():
        if graph.in_degree(n) > 1:
            list_predecesseur=[node for node in graph.predecessors(n)]
            if len(list_predecesseur) > 1:
                for i in list_predecesseur:
                    for j in list_predecesseur:
                        if i != j:
                            noeud_ancestre= nx.lowest_common_ancestor(graph,i,j)
                            if noeud_ancestre != None:
                                bubble = True
                                break
            if bubble is True:
                break
    if bubble is True:
        graph = simplify_bubbles(solve_bubble(graph,noeud_ancestre,n))
    return graph


def solve_entry_tips(graph, starting_nodes):
    graph_1 = graph.copy()
    for n in graph.nodes():
        path_list = []
        if graph.in_degree(n) > 1:
            for s in starting_nodes:
                if nx.has_path(graph,s,n) is True:
                    for p in nx.all_simple_paths(graph,s,n):
                        path_list.append(p)
            weight_avg_list=[path_average_weight(graph,path) for path in path_list]
            path_length = [len(path) for path in path_list]
            graph_1 = select_best_path(graph_1,path_list,path_length,weight_avg_list,
                delete_entry_node=True)
    return graph_1


def solve_out_tips(graph, ending_nodes):
    graph_1=graph.copy()
    for n in graph.nodes():
        path_list=[]
        if graph.out_degree(n) > 1:
            for e in ending_nodes:
                if nx.has_path(graph,n,e) is True:
                    for p in nx.all_simple_paths(graph,n,e):
                        path_list.append(p)
            weight_avg_list=[path_average_weight(graph,path) for path in path_list]
            path_length=[len(path) for path in path_list]
            graph_1=select_best_path(graph_1,path_list,path_length,weight_avg_list,
                delete_sink_node=True)
    return graph_1

def get_starting_nodes(graph):
    starting_nodes=[]
    for n in graph.nodes():
        if len(list(graph.predecessors(n))) == 0:
            starting_nodes.append(n)
    return starting_nodes

def get_sink_nodes(graph):
    sink_nodes=[]
    for n in graph.nodes():
        if len(list(graph.successors(n))) == 0:
            sink_nodes.append(n)
    return sink_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    contigs=[]
    for s in starting_nodes:
        for e in ending_nodes:
            if nx.has_path(graph,s,e) is True:
                for p in nx.all_simple_paths(graph,s,e):
                    contig=p[0]
                    for n in range(1,len(p)):
                        contig=contig+p[n][-1]
                    contigs.append((contig,len(contig)))
    return contigs

def save_contigs(contigs_list, output_file):
    f = open(output_file, "w")
    for i,contig in enumerate(contigs_list) :
        f.write(">contig_"+str(i) +" len="+str(contig[1])+"\n")
        f.write(fill(contig[0],width=80)+"\n")
    f.close()

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph"""
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
        pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)
    kmer_dict = build_kmer_dict(args.fastq_file,args.kmer_size)
    graph = build_graph(kmer_dict)

    graph = simplify_bubbles(graph)

    starting_nodes = get_starting_nodes(graph)
    graph = solve_entry_tips(graph,starting_nodes)

    ending_nodes = get_sink_nodes(graph)
    graph = solve_out_tips(graph,ending_nodes)

    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    contigs = get_contigs(graph,starting_nodes,ending_nodes)

    save_contigs(contigs,args.output_file)

    """lis=build_kmer_dict(args.fastq_file,22)
    G=build_graph(lis)
    print(G.number_of_nodes())
    start=get_starting_nodes(G)
    ending_nodes=get_sink_nodes(G)
    contigs=get_contigs(G,start,ending_nodes)
    output_file="contigs.txt"
    print(start)
    save_contigs(contigs,output_file)"""

if __name__ == '__main__':
    main()
