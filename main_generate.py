##############################################################################
# Copyright (c) 2023, Ruben Becker, Sajjad Ghobadi
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * The name of the author may not be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
# EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##############################################################################

"""Main module that executes the experiments."""

import os
import sys


import numpy as np
import networkx as nx
from numpy.random import choice, seed

import generation as gen


def create_graphs():
    """Create graph and store it in seperate files.

    """
    graphs = generator(generator_parameters)
   
    for G in graphs:
        index_node = {}
        node_index = {}
        j = 0;
        for u in sorted(G.nodes()):
            index_node[j] = u
            node_index[u] = j
            j += 1
        
        if experiment_type == 'email-Eu-core_0_0.2':
            gen.set_communities_email_Eu_core(G)
        else:
            gen.set_communities(G, comm_type)
            
        create_files_in_TIM(G)
        pass

def create_files_in_TIM(G):
     """Store the graph G in seperate files.
     
     Create seperate files for storing
         - the edges and probabilities on the edges (graph_ic.txt)
         - the number of nodes (n) and edges (m) of G (attribute.txt)
         - the community structure in G (community.txt)
     """
     
     folder = './data_set/'+ experiment_type + '/'
     path = os.path.join(folder, G.graph['graphname'])
     os.makedirs(path, exist_ok = True)
     
     graph_file = path + '/' + 'graph_ic' + '.txt'
     attribute_file = path + '/' + 'attribute' + '.txt'
     community_file = path + '/' + 'community' + '.txt'
     
     edges = sorted(G.edges(data=True), key=lambda x: x[0], reverse = False)

     with open(graph_file, 'w') as f:
         f.writelines([str(e[0]) + ' ' + str(e[1]) + ' ' + str(round(G[e[0]][e[1]]['p'], 6)) + '\n' for e in edges])
         
     with open(attribute_file, 'w') as f:
         f.writelines(['n='+ str(len(G)) + '\n', 'm='+ str(len(G.edges()))])
         
     with open(community_file, 'w') as f:
         f.writelines(str(len(G.graph['communities']))+ '\n')
         f.writelines([str(C)+ ' '+' ' .join([str(i) for i in G.graph['communities'][C]]) + '\n' for C in G.graph['communities']])

#############
# main
#############

# forbid python 2 usage
version = sys.version_info[0]
if version == 2:
    sys.exit("This script shouldn't be run by python 2 ")

# do not set seed specifically
s = None
seed(s)

# dictionary to specify different graph generators by shorter names
generators = {
    'tsang': gen.graph_tsang,
    'ba': gen.directed_barabasi_albert,
    'data_set' : gen.graph_fish,
}


print('++++++++++++++++++++++++++++++++++++++++++++++++++++')
print('++++++ Expecting experiment_type as argument. ++++++')
print('++++++++++++++++++++++++++++++++++++++++++++++++++++')

# read number of desired processes from the shell
experiment_type = sys.argv[1]
if len(sys.argv) == 3:
    number_of_processes = int(sys.argv[2])
else:
    number_of_processes = 1

rep_graph = 5


if experiment_type == 'ba-singletons-0_0.4-200':
    generator_name = 'ba'
    generator_parameters = ([200], 2, '0_0.4', rep_graph, s)
    comm_type = 'singleton' 
    
elif experiment_type == 'ba-singletons-0_0.4-100':
    generator_name = 'ba'
    generator_parameters = ([100], 2, '0_0.4', rep_graph, s)
    comm_type = 'singleton' 
    
elif experiment_type == 'ba-singletons-0_0.4-50':
    generator_name = 'ba'
    generator_parameters = ([50], 2, '0_0.4', rep_graph, s)
    comm_type = 'singleton' 
        
elif experiment_type == 'tsang-region-gender-0-0.4':
    generator_name = 'tsang'
    generator_parameters = ('0_0.4', rep_graph)
    comm_type = 'tsang_region_gender'
    
elif experiment_type == 'arena_0_0.2-random_overlap_10':
    generator_name = 'data_set'
    generator_parameters = ('0_0.2', rep_graph, ["arena.txt"])
    comm_type = 'random_overlap_10'
    
elif experiment_type == 'irvine_0_0.2_bfs_comm_10':
    generator_name = 'data_set'
    generator_parameters = ('0_0.2', rep_graph, ["irvine.txt"])
    comm_type = 'bfs_10' 
    
elif experiment_type == 'email-Eu-core_0_0.2':
    generator_name = 'data_set'
    generator_parameters = ('0_0.2', rep_graph, ["email-Eu-core.txt"])
    comm_type = 'email-Eu-core'
    
elif experiment_type == 'ca-GrQc_0_0.2-leiden':
    generator_name = 'data_set'
    generator_parameters = ('0_0.2', rep_graph, ["ca-GrQc.txt"])
    comm_type = 'leiden'  
    
elif experiment_type == 'ca-HepTh_0_0.2_random-n_10':
    generator_name = 'data_set'
    generator_parameters = ('0_0.2', rep_graph, ["ca-HepTh.txt"])
    comm_type = 'random-n_10'
    
elif experiment_type == 'facebook_combined_0_0.2_bfs_comm_2':
    generator_name = 'data_set'
    generator_parameters = ('0_0.2', rep_graph, ["facebook_combined.txt"])
    comm_type = 'bfs_2'   
      
else:
    print("Error: Unknown option.")
    assert(0)


generator = generators[generator_name]

create_graphs()
