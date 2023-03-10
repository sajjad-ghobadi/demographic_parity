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

from collections import deque

import networkx as nx
import numpy as np
from numpy import log
from numpy.random import choice


def set_to_number(G, S):
    """Compute the binary representation of the node set S.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    S : set
        set of nodes

    Returns
    -------
    int
        binary representation of S

    """
    return sum([2**i for i in range(len(G)) if list(G.nodes())[i] in S])


def number_to_set(G, i):
    """Compute the set encoded by the integer i.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    i : int
        binary representation of a set

    Returns
    -------
    list
        set of nodes encoded by i

    """
    A = list(G.nodes())
    if i > 2**(len(A)) - 1:
        print('a out of bound')
        assert(False)
    b = []
    while i != 0:
        if i & 1:
            b.append(A[0])
        i = i >> 1
        A = A[1:]
    return b


def sigma_C_p(G, p, solution_mxmn):
    """Compute expected value of sigma_C (for all C) when sampling acc to p.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    p : dict
        Dictionary with keys integers encoding sets and values probabilities.
    solution_mxmn : dict
        Dictionary with keys integers encoding sets and values probabilities.

    Returns
    -------
    dict
        for every community the average probability of its nodes being reached

    """
    if len(solution_mxmn)>0:
        comm_pro_mxmn = {S: sigma_C(G, number_to_set(G, S))
                           for S in solution_mxmn.keys()}
        
        comm_pro_mxmn_sets = {C: sum([p[0.1] * solution_mxmn[S] * comm_pro_mxmn[S][C]
                                      for S in solution_mxmn.keys()])
                              for C in G.graph['communities']}
        
        p_prime = p.copy()
        p_prime.pop(0.1)
        
        comm_probs_for_sets = {S: sigma_C(G, number_to_set(G, S))
                           for S in p_prime.keys()}
        
        comm_pro_p_sets = {C: sum([p_prime[S] * comm_probs_for_sets[S][C]
                                   for S in p_prime.keys()])
                           for C in G.graph['communities']}
        com_pro = {C: comm_pro_mxmn_sets[C]+comm_pro_p_sets[C] for C in G.graph['communities']}
    
    else:
        comm_probs_for_sets = {S: sigma_C(G, number_to_set(G, S))
                           for S in p.keys()}
        com_pro = {C: sum([p[S] * comm_probs_for_sets[S][C]
                    for S in p.keys()])
            for C in G.graph['communities']}
    return com_pro


def sigma_C_x(G, x, eps, delta):
    """Compute expected value of sigma_C (for all C) when sampling acc to x.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    x : dict
        Dictionary with keys nodes and values probabilities.
    eps : float
        eps of (eps, delta)-approximation
    delta : float
        delta of (eps, delta)-approximation

    Returns
    -------
    dict
        for every community the average probability of its nodes being reached
        when sampling sets according to x

    """
    node_probs = sigma_v_x(G, x, eps, delta)
    comm_probs = {C: 0.0 for C in G.graph['communities']}
    for v in G.nodes():
        for C in G.nodes[v]['communities']:
            comm_probs[C] += node_probs[v] / len(G.graph['communities'][C])
    return comm_probs

def sigma_p(G, p, solution_mxmn):
    """Compute the expected influence when sampling according to p.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    p : dict
        Dictionary with keys integers encoding sets and values probabilities.
    solution_mxmn : dict
        Dictionary with keys integers encoding sets and values probabilities.

    Returns
    -------
    float
        The expected influence when sampling sets according to p.

    """       
    S_sigmas_mxmn = []
    if len(solution_mxmn)>0:
        Ss_mxmn = list(solution_mxmn.keys())
        pro_mxmn = list(solution_mxmn.values())
        index_mxmn = list()
            
        for i in range(len(solution_mxmn)):
            index_mxmn.append(i)
        numm_mxmn = [choice(index_mxmn, p=pro_mxmn) for _ in range(100)] 
        sets_mxmn = []
        for j in numm_mxmn:
            sets_mxmn.append(Ss_mxmn[j])
        
        for S in sets_mxmn:
            S_sigmas_mxmn.append(sigma(G, number_to_set(G, S)))
            
        sigma_mxmn = np.mean(S_sigmas_mxmn)
            
    S_sigmas = []    
    Ss = list(p.keys())
    pro = list(p.values())
    index = list()
        
    for i in range(len(p)):
        index.append(i)
    numm = [choice(index, p=pro) for _ in range(100)] 
    sets = []
    for j in numm:
        if Ss[j] == 0.1:
            S_sigmas.append(sigma_mxmn)
        else:
            sets.append(Ss[j])
    for S in sets:
        if S == 0.0:
            S_sigmas.append(0.0)
        else:
            S_sigmas.append(sigma(G, number_to_set(G, S)))
        
    return np.mean(S_sigmas) 
 
def sigma_x(G, x):
    """Compute the expected influence when sampling according to x.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    x : dict
        Dictionary with keys nodes and values probabilities.

    Returns
    -------
    float
        The expected influence when sampling sets according to x.

    """
    N = 100
    sets = [[v for v in x.keys() if choice([0, 1], p=[1 - x[v], x[v]])]
            for i in range(N)]
    
    S_sigmas = [sigma(G, S, weight='weight') for S in sets]
    return np.mean(S_sigmas)


def sigma_v_x(G, x, eps, delta):
    """Compute (eps, delta)-apx to expected sigma_v (for all v).

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    x : dict
        Dictionary with keys nodes and values probabilities.
    eps : float
        eps of (eps, delta)-approximation
    delta : float
        delta of (eps, delta)-approximation

    Returns
    -------
    dict
        for every node the expected probability that it is reached when
        sampling sets according to x

    """
    N = int(1.0 / eps**2 * log(1.0 / delta))
    sets = [[v for v in x.keys() if choice([0, 1], p=[1 - x[v], x[v]])]
            for i in range(N)]
    S_probs = [sigma_v(G, S) for S in sets]
    return {v: 1.0 / N * sum([S_probs[i][v]
                              for i in range(len(sets))]) for v in G.nodes()}


def sigma(G, S, weight='weight'):
    """Compute the expected influence of set S.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    S : set
        The seed set of nodes.
    weight : string
        The node attribute that describes the weight of the nodes.

    Returns
    -------
    float
        The expected influence of the set S.

    """
    prob = sigma_v(G, S)
    return sum([prob[v] * G.nodes[v][weight] for v in G.nodes()])


def sigma_v(G, S):
    """Compute for all nodes the probability that they are reached by S.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    S : set
        The seed set of nodes.

    Returns
    -------
    dict
        for every node the probability that it is reached from S

    """
    node_probs = {v: 0.0 for v in G.nodes()}
    T = len(G.graph['leg'])
    for t in range(T):
        reachable_nodes = reachable(G, S, t)
        for i in reachable_nodes:
            node_probs[i] += 1.0
    return {v: node_probs[v] / T for v in G.nodes()}


def sigma_C(G, S):
    """Compute average probability of being reached in C (for all C) from S.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    S : set
        The seed set of nodes.

    Returns
    -------
    dict
        for every community the average probability of its nodes being reached
        from S

    """
    node_probs = sigma_v(G, S)
    comm_probs = {C: 0.0 for C in G.graph['communities']}
    for v in G.nodes():
        for C in G.nodes[v]['communities']:
            comm_probs[C] += node_probs[v] / len(G.graph['communities'][C])
    return comm_probs


def live_out_neighbors(G, v, outcome):
    """Return the live out neighbors of v in outcome.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    v : int
        Node whose neighbors are seeked.
    outcome : list
        list of edges that are live

    Returns
    -------
    list
        List of out neighbors that are incident to v via a live edge

    """
    return [w for w in list(G.successors(v)) if (v, w) in outcome]


def reachable(G, seed_set, t):
    """Compute a set of nodes reachable from seed_set in outcome t.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    seed_set : list
        List of seed nodes.
    t : int
        Integer identifying which outcome is asked.

    Returns
    -------
    set
        Set of nodes reachable from seed_set in outcome t.

    """
    reachable = set(seed_set)
    if 'reachable_nodes' in G.graph:
        for u in seed_set:
            reachable |= G.graph['reachable_nodes'][t][u]
    else:
        queue = deque([(seed, live_out_neighbors(G, seed, G.graph['leg'][t]))
                       for seed in seed_set])
        while queue:
            parent, children = queue.popleft()
            for child in children:
                if child not in reachable:
                    reachable.add(child)
                    queue.append((child,
                                  live_out_neighbors(G,
                                                     child,
                                                     G.graph['leg'][t])))
    return reachable


def uniform(G):
    """Return the uniform node based solution, i.e., k/n for every node.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.

    Returns
    -------
    dict
        For every node in G, the uniform probability k/n.

    """
    return {i: G.graph['k'] / len(G) for i in G.nodes()}

def gen_leg(G, T=None):
    """Generate T live-edge graphs according to the IC model.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.

    Returns
    -------
    pair
        of live-edge-graphs as list of edges
        and as networkx.Graphs (for the implementation of Tsang et al.)

    """
    if T is None:
        T = G.graph['T']
    legs, leg_tsang = [], []
    for ell in range(T):
        H = nx.Graph()
        H.add_nodes_from(sorted(G.nodes()))
        legs.append(set())
        for u, v in sorted(G.edges()):
            if choice([0, 1], p=[1 - G[u][v]['p'], G[u][v]['p']]):
                legs[ell].add((u, v))
                H.add_edge(u, v)
        leg_tsang.append(H)
    return legs, leg_tsang


def compute_marginals_im(G, V_w_S, reachable_from_s, weight='weight'):
    """Compute the marginal increment to sigma of all nodes in V_w_S.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    V_w_S : set
        Set of nodes that are not in S.
    reachable_from_s : set
        Set of nodes that are reachable from S.
    weight : string
        The node attribute that describes the weight of the nodes.

    Returns
    -------
    dict
        For every node in V_w_S the marginal increment it can give
        when added to S.

    """
    return {v: 1 / G.graph['T'] * sum([sum([G.nodes[w][weight]
                                            for w in
                                            G.graph['reachable_nodes'][t][v]
                                            if w not in reachable_from_s[t]])
                                       for t in range(G.graph['T'])])
            for v in V_w_S}


def greedy_influence_max(G, weight='weight'):
    """Greedy algorithm for influence maximization.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    weight : string
        The node attribute that describes the weight of the nodes.

    Returns
    -------
    pair
        of computed set S of nodes and the objective value achieved by S.

    """
    assert(G.graph['k'] <= len(G))

    if 'reachable_nodes' not in G.graph:
        G.graph['reachable_nodes'] = {t: {v: reachable(G,
                                                       set([v]),
                                                       t)
                                          for v in G.nodes()}
                                      for t in range(G.graph['T'])}

    S = set()
    V_w_S = set(list(G.nodes()))

    reachable_from_s = {t: set([]) for t in range(G.graph['T'])}

    increment_sum = 0

    while len(S) < G.graph['k']:
        marginals = compute_marginals_im(G, V_w_S,
                                         reachable_from_s, weight=weight)

        max_v = max(marginals, key=lambda v: marginals[v])

        increment_sum += marginals[max_v]
        S.add(max_v)
        V_w_S.remove(max_v)

        reachable_from_s = {t: G.graph['reachable_nodes'][t][max_v]
                            | reachable_from_s[t] for t in range(G.graph['T'])}

    return S, sigma(G, S, weight=weight)


def compute_s_i_probs(G, V_w_S, reachable_from_s):
    """Compute sigma_C(S u i) for all i and C.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    V_w_S : set
        Set of nodes that are not in S.
    reachable_from_s : set
        Set of nodes that are reachable from S.

    Returns
    -------
    pair
        of
        dict r_s_i: set of reachable nodes from S u {i} for each t in [T]
        dict: for each i in V_w_S the minimum average probability of nodes to
        be reached from S u {i} of nodes in communities C

    """
    r_s_i = {i: {
        t: reachable_from_s[t] | G.graph['reachable_nodes'][t][i]
        for t in range(G.graph['T'])}
        for i in V_w_S}  
    
    probs_v = {i: {v: 1 / G.graph['T'] * sum([int(v in x)
                                           for x in r_s_i[i].values()])
                   for v in G.nodes()}
               for i in V_w_S}
    
    comms = G.graph['communities']
    
    return r_s_i, {i:
                    min([1 / len(comms[C]) * sum([probs_v[i][v]
                                                  for v in comms[C]])
                        for C in comms.keys()])
                    for i in V_w_S}


def grdy_maximin(G, weight='weight'):
    """Implement the greedy maximin algorithm.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    weight : string
        The node attribute that describes the weight of the nodes.

    Returns
    -------
    set
        Computed seed set.

    """

    if 'reachable_nodes' not in G.graph:
        G.graph['reachable_nodes'] = {t: {v: reachable(G,
                                                       set([v]),
                                                       t)
                                          for v in G.nodes()}
                                      for t in range(G.graph['T'])}
        
    
    V_w_S = set(list(G.nodes()))
    max_degree = max(G.out_degree(weight='p'), key=lambda x: x[1])[0]
    S = set([max_degree])

    V_w_S.remove(max_degree)

    reachable_from_s = {t: reachable(G, set([max_degree]), t)
                            for t in range(G.graph['T'])}
 
 
    while len(S) < G.graph['k']:
        
        r_s_i, min_probs_s_i = compute_s_i_probs(
            G, V_w_S, reachable_from_s)
                
        max_v = max(min_probs_s_i, key=lambda i: (min_probs_s_i[i], i))
    
        S.add(max_v)
        V_w_S.remove(max_v)
    
        reachable_from_s = r_s_i[max_v]

    return S
