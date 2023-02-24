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

"""Module proxy for the cpp implementation."""

import os
import uuid
import influence_max as im
import operator
import random
import time
from decimal import Decimal
  


def grdy_prop(G, graph_file, k, num):
    """Call execute_deterministics for executing the grdy_prop from the cpp implementation.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    graph_file : The path of the folder containing graph files.
    k : int
        Budget (size of the seed set)    
    num : int
        Number specifying the name of output file using cpp implementation.

    """
    return execute_deterministics(G, graph_file, k, num, 'greedy') 


def grdy_im(G, graph_file, k, num):
    """Call execute_deterministics for executing the grdy_im from the cpp implementation.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    graph_file : The path of the folder containing graph files.
    k : int
        Budget (size of the ssed set)    
    num : int
        Number specifying the name of output file using cpp implementation.

    """
    return execute_deterministics(G, graph_file, k, num, 'tim') 


def execute_deterministics(G, graph_file, k, num, setting): 
    """Execute deterministic algorithms from the cpp implementation.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    graph_file : The path of the folder containing graph files.
    k : int
        Budget (size of the ssed set)    
    num : int
        Number specifying the name of output file using cpp implementation.
    setting : str
        String specifying which function should be executed.
        
    Returns
    -------
    S : list of nodes
        Computed solution.
    running_time : The execution time of the function (according to setting) to return the computed solution.
    
    """
    log_file = str(random.random()) + str(uuid.uuid1()) + str(num) + '.txt'
    
    if setting == 'tim':
        path_make = 'cd' + ' ' + './TIM_final' + '&& make'
        tim_command = 'cd' + ' ' + './TIM_final' + '&& ./node_set_based -model IC -dataset' + ' ' +'.'+ graph_file + ' ' + '-epsilon 0.1 -k' + ' ' + str(k) + ' ' + '-setting' + ' '+ 'tim' + ' ' + '-file_name' + ' ' + str(log_file)
        file = './TIM_final/' + log_file
    else:
        path_make = 'cd' + ' ' + './TIM_final' + '&& make'
        tim_command = 'cd' + ' ' + './TIM_final' + '&& ./node_set_based -model IC -dataset' + ' ' +'.'+ graph_file + ' ' + '-epsilon 0.1 -k' + ' ' + str(k) + ' ' + '-setting' + ' '+ 'greedy' + ' ' + '-file_name' + ' ' + str(log_file)
        file = './TIM_final/' + log_file
     
    process_make = os.popen(path_make)
    output_make = process_make.read()
    process_make.close()
    print(output_make)

    process = os.popen(tim_command)
    output = process.read()
    process.close()
    print(output) 
    S = set()
    with open(file) as f:
        line1 = f.readline()
        running_time = float(line1)
        
        line2 = f.readline()
        nodes = map(str, line2.split())
        for v in nodes:
            S.add(int(v))

    os.remove(file)
    return S, running_time


def mult_weight(G, graph_file, k, num):
    """Execute multiplicative weight routine for the set-based problem of Becker et al (2022) from cpp implementation.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    graph_file : The path of the folder containing graph files.
    k : int
        Budget (size of the ssed set)    
    num : int
        Number specifying the name of output file using cpp implementation.
        
    Returns
    -------
    p : dict
        dictionary with keys integers encoding sets and values probabilities.
    running_time : The execution time of mult_weight to return the computed solution

    """
    path_make = 'cd' + ' ' + './TIM_final' + '&& make'
    process_make = os.popen(path_make)
    output_make = process_make.read()
    process_make.close()
    print(output_make)
    
    log_file = str(random.random()) + str(uuid.uuid1()) + str(num) + '.txt'
    tim_command = 'cd' + ' ' + './TIM_final' + '&& ./node_set_based -model IC -dataset' + ' ' +'.' + graph_file + ' ' + '-epsilon 0.1 -k' + ' ' + str(k) + ' ' + '-setting' + ' '+ 'set' + ' ' + '-file_name' + ' ' + str(log_file)
    file = './TIM_final/'+ log_file
    process = os.popen(tim_command)
    output = process.read()
    process.close()
    print(output)   
        
    p = {}
    with open(file) as f:
        line1 = f.readline()
        running_time = float(line1)
            
    with open(file) as f:
        lines = f.readlines()[1:]
        for line in lines:
            *sets , pro = map(float, line.split())
            if len(sets) == 0.0:
                p[0.0] = pro
            else:
                p[im.set_to_number(G, sets)] = pro
    os.remove(file)
    sum_p = sum(p.values())
    if sum_p != 1:
        p = dict(sorted(p.items(), key=lambda x: x[1]))
        if sum_p > 1:
            x = sum_p - 1
            for S_pro in p:
                if p[S_pro] >= x:
                      p[S_pro] -= x
                      break
        if sum_p < 1:
            x = 1 - sum_p
            max_S = max(p.items(), key=operator.itemgetter(1))[0]
            p[max_S] += x;                       
    return p, running_time


def grdy_grp_lp(G, graph_file, k, num, folder):
    """Call execute_cpp for executing the grdy_grp_lp in the paper from the cpp implementation.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    graph_file : The path of the folder containing graph files.
    k : int
        Budget (size of the ssed set)    
    num : int
        Number specifying the name of output file using cpp implementation.
    folder : str
        The path of the folder for writing the computed solutions in seperate files.

    """
    execute_cpp(G, graph_file, k, num, folder, 'set_dp2')
     

def ind_lp(G, graph_file, k, num, folder):
    """Call execute_cpp for executing the ind_lp in the paper from the cpp implementation.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    graph_file : The path of the folder containing graph files.
    k : int
        Budget (size of the ssed set)    
    num : int
        Number specifying the name of output file using cpp implementation.
    folder : str
        the path of the folder for writing the computed solutions in seperate files.

    """
    execute_cpp(G, graph_file, k, num, folder, 'node_dp')


def maximin_lp(G, graph_file, k, num, folder):
    """Call execute_cpp for executing the maximin_lp in the paper from the cpp implementation.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    graph_file : The path of the folder containing graph files.
    k : int
        Budget (size of the ssed set)    
    num : int
        Number specifying the name of output file using cpp implementation.
    folder : str
       The path of the folder for writing the computed solutions in seperate files.

    """
    execute_cpp(G, graph_file, k, num, folder, 'set_dp3')
        
    
def execute_cpp(G, graph_file, k, num, folder, setting):
    """Execute probabilistic algorithms in the paper from the cpp implementation.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    graph_file : The path of the folder containing graph files.
    k : int
        Budget (size of the ssed set)    
    num : int
        Number specifying the name of output file using cpp implementation.
    folder : str
        The path of the folder for writing the computed solutions in seperate files.
    setting : str
            string specifying which function should be executed.

    """
    path_make = 'cd' + ' ' + './TIM_final' + '&& make'
    process_make = os.popen(path_make)
    output_make = process_make.read()
    process_make.close()
    print(output_make)
    
    log_file = str(time.time()) + str(uuid.uuid1()) + '_' + str(num)
    tim_command = 'cd' + ' ' + './TIM_final' + '&& ./node_set_based -model IC -dataset' + ' ' +'.' + graph_file + ' ' + '-epsilon 0.1 -k' + ' ' + str(k) + ' ' + '-setting' + ' '+ setting + ' ' + '-file_name' + ' ' + str(log_file) + ' ' + '-folder' + ' ' + str(folder)
    
    process = os.popen(tim_command)
    output = process.read()
    process.close()
    print(output) 
    
    
def read_file_dp3(G, file):
    """Read the solution written in file after executing maximin_lp.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    file : file containing the computed solution.
    
    Returns
    -------
    p : dict
        dictionary with keys integers encoding sets and values probabilities.
    p_maximin : dict
        dictionary with keys integers encoding sets and values probabilities (the solution returned by the algorithm of Becker et al 2022).
    omega : str
        violation in fairness
    epss2 : double
        violation in fairness
    running_time : The execution time of maximin_lp to return the computed solution

    """
    p = {}
    with open(file) as f:
        line1 = f.readline()
        running_time = float(line1)
        
        line2 = f.readline()
        omega = str.strip(line2)
        
        line3 = f.readline()
        epss2 = str.strip(line3)
            
    with open(file) as f:
        lines = f.readlines()[3:]
        counter = 0
        for line in lines:
            counter += 1  
            line2 = line.rstrip('\n')
            if line2 == 'maximin':
                break
            if line2 != 'maximin':
                *sets , pro = map(float, line.split())
                if len(sets) == 0.0:
                    p[0.0] = pro
                else:
                    p[im.set_to_number(G, sets)] = pro
        
    with open(file) as f:
        a = 0
        p_maximin = {}
        lines_mxmn = f.readlines()[3 + counter:]
        for line in lines_mxmn:
            a += 1
            if a == 1:
                pro_maximin = float(line)
            if a>1:
                *sets , pro = map(float, line.split())
                p_maximin[im.set_to_number(G, sets)] = pro
            
    os.remove(file)    
    
    if len(p_maximin) > 0:        
        sum_p_mxmn = sum(p_maximin.values())
        if sum_p_mxmn != 1:
            p_maximin = dict(sorted(p_maximin.items(), key=lambda x: x[1]))
            if sum_p_mxmn > 1:
                x = sum_p_mxmn - 1
                for S_pro in p_maximin:
                    if p_maximin[S_pro] >= x:
                          p_maximin[S_pro] -= x
                          break
            if sum_p_mxmn < 1:
                x = 1 - sum_p_mxmn
                max_S = max(p_maximin.items(), key=operator.itemgetter(1))[0]
                p_maximin[max_S] += x; 
    
        p[0.1] = pro_maximin  
    
    
    sum_p = sum(p.values())
    if sum_p != 1:
        p = dict(sorted(p.items(), key=lambda x: x[1]))
        if sum_p > 1:
            x = sum_p - 1
            for S_pro in p:
                if p[S_pro] >= x:
                      p[S_pro] -= x
                      break
        if sum_p < 1:
            x = 1 - sum_p
            max_S = max(p.items(), key=operator.itemgetter(1))[0]
            p[max_S] += x;  
                      
    return p, p_maximin, omega, epss2, running_time

        
def read_file_dp2(G, file):
    """Read the solution written in file after executing grdy_grp_lp.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    file : file containing the computed solution.
    
    Returns
    -------
    p : dict
        dictionary with keys integers encoding sets and values probabilities.
    omega : str
        violation in fairness
    epss2 : double
        violation in fairness
    running_time : The execution time of grdy_grp_lp to return the computed solution

    """
    p = {}
    with open(file) as f:
        line1 = f.readline()
        running_time = float(line1)
        
        line2 = f.readline()
        omega = str.strip(line2)
        
        line3 = f.readline()
        epss2 = float(line3)
            
    with open(file) as f:
        lines = f.readlines()[3:]
        for line in lines:
            *sets , pro = map(float, line.split())
            if len(sets) == 0.0:
                p[0.0] = pro
            else:
                p[im.set_to_number(G, sets)] = pro
    os.remove(file)
    sum_p = sum(p.values())
    if sum_p != 1:
        p = dict(sorted(p.items(), key=lambda x: x[1]))
        if sum_p > 1:
            x = sum_p - 1
            for S_pro in p:
                if p[S_pro] >= x:
                      p[S_pro] -= x
                      break
        if sum_p < 1:
            x = 1 - sum_p
            max_S = max(p.items(), key=operator.itemgetter(1))[0]
            p[max_S] += x;                       
    return p, omega, epss2, running_time


def read_node_dp(G, file):
    """Read the solution written in file after executing ind_lp.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    file : file containing the computed solution.
    
    Returns
    -------
    x : dict
        dictionary with keys nodes and values probabilities.
    omega : str
        violation in fairness
    epss2 : double
        violation in fairness
    running_time : The execution time of ind_lp to return the computed solution

    """
    x = {}
    with open(file) as f:
        line1 = f.readline()
        running_time = float(line1)
        
        line2 = f.readline()
        omega = str.strip(line2)
        
        line3 = f.readline()
        epss2 = float(line3)
        
    with open(file) as f:
        lines = f.readlines()[3:]
        for line in lines:
            node , pro = map(float, line.split())
            x[node] = pro
    os.remove(file)
    return x, omega, epss2, running_time 
 

