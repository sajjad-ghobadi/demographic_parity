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
import threading
import time
import uuid
import random

import networkx as nx
import numpy as np
from numpy.random import choice, seed

import influence_max as im
import cpp_proxy as cpp
import print_functions as pf
import generation as gen
import gurobi_influence_max as gim
import maximin_fish as mf
import moso as tm
from tsang.algorithms import rounding


def sample_sets(G, vector, times, type):
    """Sample times many sets from vector depending on type.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    vector : depends on type
        If type is deterministic, vector is a set.
        If type is node_based or swap_rounding, vector is a dictionary with
            keys being the nodes and values the probability of the respective
            node.
        If type is set_based, vector is a dictionary with keys being ints
            representing the sets (binary) and values are the probability
            of the respective set.
    times : int
        How many times to sample a set.
    type : str
        The type that specifies how to sample the sets, this depends on whether
        this is a deterministic or probabilistic problem.

    Returns
    -------
    list of lists
        The list of sampled sets.

    """
    if type == 'deterministic':
        sets = [vector for _ in range(times)]
    elif type == 'node_based':
        sets = [[v for v in vector.keys() if choice(
            [0, 1], p=[1 - vector[v], vector[v]])]
            for _ in range(times)]
    elif type == 'swap_rounding':
        x_items_list = sorted(vector.items())
        x = np.array([x_items[1] for x_items in x_items_list])
        rounded_xs = [rounding(x) for _ in range(times)]
        sets = [[v for v in G.nodes() if rounded_xs[i][v]]
                for i in range(times)]
    elif type == 'set_based':
        Ss = list(vector.keys())
        pro = list(vector.values())
        index = list()
        
        for i in range(len(vector)):
            index.append(i)
        numm = [choice(index, p=pro) for _ in range(times)] 
        sets = []
        for j in numm:
            sets.append(Ss[j])                                
    else:
        print("Error: Unknown option.", type)
        assert(False)
    return sets


def comp_ex_post(G, solution, fct_name, solution_mxmn):
    """Compute ex_post values by sampling one set.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    solution : depends on fct_name
        The computed solution.
    fct_name : str
        The name of the function that computed solution.
    solution_mxmn : dict
        The solution of the algorithm of Becker et al (2022) returned from maximin_lp.
        
    Returns
    -------
    float
        The ex_post fairness value and spread obtained for one sampled set.

    """ 
    sets = sample_sets(G, solution, 1,
                        ex_post_sampling_types[fct_name])
    min_probs_sets = []
    max_probs_sets = []
    
    if fct_name == 'grdy_grp_lp' or fct_name == 'mult_weight' or fct_name == 'maximin_lp':
        for S in sets:       
            if S == 0.1:
                sets.clear()
                Ss = list(solution_mxmn.keys())
                pro = list(solution_mxmn.values())
                index = list()
                
                for i in range(len(solution_mxmn)):
                    index.append(i)
                numm = [choice(index, p=pro) for _ in range(1)] 
                sets = []
                for j in numm:
                    sets.append(Ss[j])
                
        sigma_sets_ep = [im.sigma(G, im.number_to_set(G, S), weight='weight') for S in sets]
        for S in sets:
            pro_C = im.sigma_C(G, im.number_to_set(G, S)).values()
            min_probs_sets.append(min(pro_C))
            max_probs_sets.append(max(pro_C))
    else:
        sigma_sets_ep = [im.sigma(G, S, weight='weight') for S in sets]
        for S in sets:
            pro_C = im.sigma_C(G, S).values()
            min_probs_sets.append(min(pro_C))
            max_probs_sets.append(max(pro_C))
    sigma_ep = np.mean(sigma_sets_ep)    
    ex_post_min = np.mean(min_probs_sets)
    ex_post_max = np.mean(max_probs_sets)
    
    return ex_post_min, ex_post_max, sigma_ep


def comp_ex_ante(G, solution, fct_name, solution_mxmn):
    """Compute ex_ante values.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    solution : depends on fct_name
        The computed solution.
    fct_name : str
        The name of the function that computed solution.
    solution_mxmn : dict
        The solution of algorithm of Becker et al (2022) returned from maximin_lp.
    
    Returns
    -------
    float
        The ex_ante fairness value and the expected spread.

    """
    if fct_name in deterministic_algorithms + ['grdy_im'] + ['grdy_prop'] + ['milp']:
        comm_probs = im.sigma_C(G, solution)
        return min(comm_probs.values()), max(comm_probs.values()), 0.0
    else:
        assert(fct_name in probabilistic_algorithms + cpp_algos)
        if fct_name in ['mult_weight', 'grdy_grp_lp','maximin_lp']:
            comm_probs = im.sigma_C_p(G, solution, solution_mxmn)  
            sigma = im.sigma_p(G, solution, solution_mxmn)
        elif fct_name in ['ind_lp', 'moso', 'uniform']:
            comm_probs = im.sigma_C_x(G, solution, 0.1, 0.1)   
            sigma = im.sigma_x(G, solution)
        else:
            print("Error: Unknown option:", fct_name)
            assert(False)
        return min(comm_probs.values()) , max(comm_probs.values()), sigma


def read_graph(graph_file):
    """Read graph from files.

    Parameters
    ----------
    graph_file : The path of the folder containing graph files

    Returns
    -------
    G : networkx.DiGraph
        The underlying considered instance.
    """
    G = nx.read_edgelist(graph_file + '/graph_ic.txt', nodetype = int,
                         data=(('p',float),), create_using = nx.DiGraph())
      
    gen.init_communities(G)    
    G.graph['nodes_in_comms'] = set()    
    
    comm_file = open(graph_file + '/community.txt','r')
    lines = comm_file.readlines()[1:]    
    for line in lines:
        comm , *nodes = map(str, line.split())
        G.graph['communities'][comm] = [int(v) for v in nodes]
        
        for u in nodes:
            G.nodes[int(u)]['communities'].append(comm)
            G.graph['nodes_in_comms'].add(int(u))        
    
    comm_file.close()    
    gen.set_unit_node_weights(G)
             
    return G


def execute(function, graph_file, k, rep_per_graph_i, num): 
    """Execute function on G (using graph_file) and writes results to out_file (global).

    Parameters
    ----------
    function : function that takes a networkx graph and returns a solution.
    graph_file : The path of the folder containing graph files.
    k : int
        Budget (size of the seed set)
    rep_per_graph_i : int
        The number showing repetition per graph.
    
    num : int
        Number specifying the name of output file (used for cpp implementations).

    """
    fct_name = function.__name__
        
    G = read_graph(graph_file)
    G.graph['T'] = no_le_graphs
    G.graph['k'] = k
            
    G.graph['leg'], G.graph['leg_tsang'] = im.gen_leg(G)  
    start = time.time()
    if fct_name in python_algos:
        solution = function(G)
        ex_time = time.time() - start
        final_evaluation(G, solution, fct_name, graph_file, k, ex_time, rep_per_graph_i, 0.0, 0.0, [])
    else:
        if fct_name in ['ind_lp', 'grdy_grp_lp', 'maximin_lp']:
            folder_name = str(num) + '_' + str(time.time()) + experiment_type + str(uuid.uuid1())
            if fct_name == 'ind_lp':
                folder = './TIM_final/results_node_dp' + folder_name
                folder_cpp = './results_node_dp' + folder_name + '/'
                os.mkdir(folder)
                function(G, graph_file, k, num, folder_cpp)
                file_list = sorted([f for f in os.listdir(folder+'/') if f.endswith(".txt")])
                for file in file_list:
                    solution, omega, epss2, ex_time = cpp.read_node_dp(G, folder+'/'+file)
                    final_evaluation(G, solution, fct_name, graph_file, k, ex_time, rep_per_graph_i, epss2, omega, [])
                os.rmdir(folder)
                
            elif fct_name == 'grdy_grp_lp':
                folder = './TIM_final/results_dp2' + folder_name
                folder_cpp = './results_dp2' + folder_name + '/'  
                
                os.mkdir(folder)
                function(G, graph_file, k, num, folder_cpp)
                file_list = sorted([f for f in os.listdir(folder+'/') if f.endswith(".txt")])
                for file in file_list:
                    solution, omega, epss2, ex_time = cpp.read_file_dp2(G, folder+'/'+file)        
                    final_evaluation(G, solution, fct_name, graph_file, k, ex_time, rep_per_graph_i, epss2, omega, [])
                os.rmdir(folder)
                    
            else:
                fct_name == 'maximin_lp'
                folder = './TIM_final/results_dp3' + folder_name
                folder_cpp = './results_dp3' + folder_name + '/'            
                os.mkdir(folder)
        
                function(G, graph_file, k, num, folder_cpp)
                file_list = sorted([f for f in os.listdir(folder+'/') if f.endswith(".txt")])
                for file in file_list:
                    solution, solution_mxmn, omega, epss2, ex_time = cpp.read_file_dp3(G, folder+'/'+file)        
                    final_evaluation(G, solution, fct_name, graph_file, k, ex_time, rep_per_graph_i, epss2, omega, solution_mxmn)
                os.rmdir(folder)
                                    
        else:
            solution, ex_time = function(G, graph_file, k, num)
            final_evaluation(G, solution, fct_name, graph_file, k, ex_time, rep_per_graph_i, 0.0, 0.0, [])


def final_evaluation(G, solution, fct_name, graph_file, k, ex_time, rep_per_graph_i, epss, omega, solution_mxmn):
    """Compute objective values using the solution and write results to out_file (global).

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    solution : depends on fct_name
        The computed solution.
    fct_name : str
        The name of the function that computed solution.
    graph_file : The path of the folder containing graph files.
    k : int
        Budget (size of the ssed set)
    ex_time : The execution time of function to return solution.
    rep_per_graph_i : int
        The number showing repetition per graph.
    epss : float
        violation in fairness
    omega : str
        violation in fairness
    solution_mxmn : dict
        The solution of the algorithm of Becker et al (2022) returned from maximin_lp.

    """
    
    G.graph['leg'], G.graph['leg_tsang'] = im.gen_leg(G) 
    n, m = len(G), len(G.edges())
    ex_ante_min, ex_ante_max, sigma_ea = comp_ex_ante(G, solution, fct_name, solution_mxmn) 
    ex_post_min, ex_post_max, sigma_ep = comp_ex_post(G, solution, fct_name, solution_mxmn)
    x = graph_file.split("/")
    graph_file = x[3]
    
    if fct_name in ['grdy_grp_lp', 'maximin_lp', 'ind_lp']:
        fct_name = fct_name + '_' + str(omega)
    if fct_name in deterministic_algorithms + ['grdy_im'] + ['grdy_prop'] + ['milp']:
        sigma_ea = sigma_ep
    
    res = (graph_file,
           fct_name,
           n,
           m,
           k,
           ex_time,
           rep_per_graph_i,
           omega,
           epss,
           sigma_ep,
           sigma_ea,
           ex_post_min,
           ex_post_max,
           ex_ante_min,
           ex_ante_max)
    
    with open(out_file, 'a') as f:
        f.write('\t\t\t'.join('%s' % x for x in res) + '\n')

    print("\n")
    print("function: ".ljust(30, ' '), fct_name)
    # print("maximizing solution: ".ljust(30, ' '))
    # if fct_name == 'mult_weight' or fct_name == 'grdy_grp_lp' or fct_name == 'maximin_lp':
    #     pf.myprint([(im.number_to_set(G, s), solution[s])
    #                 for s in solution.keys() if s!=0])
    # else:
    #     pf.myprint(solution)
    print("ex_ante_min: ".ljust(30, ' '), ex_ante_min)
    print("ex_ante_max: ".ljust(30, ' '), ex_ante_max)
    print("ex_post_min: ".ljust(30, ' '), ex_post_min)
    print("ex_post_max: ".ljust(30, ' '), ex_post_max)
    print("sigma: ".ljust(30, ' '), sigma_ea)
    print("running time: ".ljust(30, ' '), ex_time)
    if fct_name not in cpp_algos:
        print("number of live edge graphs: ".ljust(30, ' '), G.graph['T'])
    print("\n")


def generate_executables(
        functions,
        instance_folder,
        k_range,
        rep_per_graph,
        no_le_graphs):     
    """Generate executables by generating graphs and combining them.

    Parameters
    ----------
    functions : list
        List of all the functions that are supposed to be executed.
    instance_folder : the folder of the instance to be executed
    k_range : list
        Range of k-values to be tested.
    rep_per_graph : int
        Number of times the execution should be repeated per graph.
    no_le_graphs : int
        Number of live-edge graphs to be generated and used in all following
        computations.

    Returns
    -------
    list
        List of executables, i.e., (function, graph)-pairs.

    """
    num = 0     
    executables = []
    instance_names = sorted([f for f in os.listdir('./data_set/' + instance_folder) 
                            if not f.startswith(".")])
    for graph_folder_name in instance_names:   
        for rep_per_graph_i in range(rep_per_graph):
            for k in k_range:
                for function in functions:
                    num += 1
                    executables.append((function,
                                        './data_set/' + instance_folder + '/' + graph_folder_name,
                                            k, rep_per_graph_i, num)) 
            print('collected', len(executables), 'executables')
    return executables


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

# the following lists are used to specify the way in which sets are sampled
# for the respective algorithms

deterministic_algorithms = [
    'myopic',
    'naive_myopic',
    'grdy_maximin',
    'milp']

probabilistic_algorithms = [
    'uniform',
    'moso']

python_algos = deterministic_algorithms + probabilistic_algorithms
cpp_algos = [
    'ind_lp',
    'grdy_grp_lp',
    'maximin_lp',
    'mult_weight']


sampling_types = [
    'deterministic',
    'node_based',
    'set_based',
    'swap_rounding'
]

ex_ante_sampling_types = {
    'mult_weight': 'set_based',
    'grdy_grp_lp': 'set_based',
    'maximin_lp': 'set_based',
    'ind_lp': 'node_based',
    'uniform': 'node_based',
    'moso': 'node_based'}
for alg in deterministic_algorithms + ['grdy_im'] + ['grdy_prop'] + ['milp']:
    ex_ante_sampling_types[alg] = 'deterministic'

ex_post_sampling_types = ex_ante_sampling_types
ex_post_sampling_types['moso'] = 'swap_rounding'

print('++++++++++++++++++++++++++++++++++++++++++++++++++++')
print('++++++ Expecting experiment_type as argument. ++++++')
print('++++++++++++++++++++++++++++++++++++++++++++++++++++')


# read number of desired processes from the shell
experiment_type = sys.argv[1]
if len(sys.argv) == 3:
    number_of_processes = int(sys.argv[2])
else:
    number_of_processes = 1

# default values for experiments     
functions = [cpp.grdy_grp_lp,
             cpp.maximin_lp,
             cpp.mult_weight,
             cpp.grdy_im,
             gim.milp,
             cpp.grdy_prop,
             im.uniform,
             tm.moso,
             mf.myopic,
             im.grdy_maximin,
             cpp.ind_lp
             ]
no_le_graphs = 100
rep_per_graph = 10

# specify experiment dependent parameters
if experiment_type == 'ba-singletons-0_0.4-200':
    instance_folder = 'ba-singletons-0_0.4-200'
    k_range = [25]
            
elif experiment_type == 'ba-singletons-0_0.4-100':
    instance_folder = 'ba-singletons-0_0.4-100'
    k_range = [25]
    
elif experiment_type == 'ba-singletons-0_0.4-50':
    instance_folder = 'ba-singletons-0_0.4-50'
    k_range = [25]
                
elif experiment_type == 'tsang-region-gender-0-0.4':
    instance_folder = 'tsang-region-gender-0-0.4'
    k_range = [25]     
        
elif experiment_type == 'arena_0_0.2-random_overlap_10':           
    instance_folder = 'arena_0_0.2-random_overlap_10'
    k_range = [100]
               
elif experiment_type == 'irvine_0_0.2_bfs_comm_10':           
    instance_folder = 'irvine_0_0.2_bfs_comm_10'
    k_range = [50]
    
elif experiment_type == 'email-Eu-core_0_0.2':           
    instance_folder = 'email-Eu-core_0_0.2'
    k_range = [100]
       
elif experiment_type == 'ca-GrQc_0_0.2-leiden':           
    instance_folder = 'ca-GrQc_0_0.2-leiden'
    k_range = [100]
            
elif experiment_type == 'ca-HepTh_0_0.2_random-n_10':           
    instance_folder = 'ca-HepTh_0_0.2_random-n_10'
    k_range = [100]
    
elif experiment_type == 'facebook_combined_0_0.2_bfs_comm_2':           
    instance_folder = 'facebook_combined_0_0.2_bfs_comm_2'
    k_range = [50]
      
else:
    print("Error: Unknown option.")
    assert(0)

# create output file with header
folder = './results_dp/'
out_file = folder + experiment_type + '.txt'
if os.path.exists(out_file):
    out_file = out_file[:-4] + '_' + str(int(time.time())) + '.txt'
print('Output is written to:', out_file, '\n')
header = [
    'graphname',
    'algorithm',
    'n',
    'm',
    'k',
    'running_time',
    'rep_per_graph',
    'omega',
    'epss',
    'spread',
    'exp_spread',
    'ex_post_min',
    'ex_post_max',
    'ex_ante_min',
    'ex_ante_max',]
with open(out_file, 'a') as f:
    f.write('\t\t\t'.join('%s' % x for x in header) + '\n')

#generate the various experiments
executables = generate_executables(
    functions,
    instance_folder,
    k_range,
    rep_per_graph,
    no_le_graphs)

# run experiments in parallel (if number_of_processes > 1)
thread_list = []
for executable in executables:
    thread_list.insert(
        0,
        threading.Thread(
            target=execute,
            args=(executable[0], executable[1], executable[2],
                  executable[3], executable[4]))) 

while thread_list:
    if threading.active_count() < number_of_processes + 1:
        thread_list.pop().start()
        