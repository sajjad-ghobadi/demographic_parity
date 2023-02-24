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

import gurobipy as gp
import influence_max as im


def milp(G):
    """ Implement the MILP-based algorithm of Farnadi et al 2020.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
        
    Returns
    -------
    set of nodes

    """
    m = gp.Model()
    m.setParam('OutputFlag', False) 
    x, y = {}, {}
    for i in G.nodes():
        x[i] = m.addVar(lb=0, ub=1, vtype=gp.GRB.BINARY, name="x%d" % i)

    for t in range(len(G.graph['leg'])):
        for j in G.nodes():
            y[(t, j)] = m.addVar(lb=0, ub=1,
                                 vtype=gp.GRB.CONTINUOUS, name="y%d,%d" % (t, j))
    
    m.update()
    m.addConstr(sum(x[i] for i in G.nodes()) <= G.graph['k'])
    for j in G.nodes():
        for t in range(len(G.graph['leg'])):
            m.addConstr(sum(x[i] for i in G.nodes() if j in im.reachable(
                G, [i], t)) >= y[(t, j)])
    
    T = len(G.graph['leg'])
    for C in G.graph['communities']:
        m.addConstr(sum(1/T * sum(y[(t, v)] for v in G.graph['communities'][C])
                        - 1/T * len(G.graph['communities'][C])/len(G.nodes()) * sum(y[(t, u)] 
                                                                                    for u in G.nodes()) 
                        for t in range(T)) <= 0.1)
        
        m.addConstr(sum(1/T * sum(y[(t, v)] for v in G.graph['communities'][C])
                        - 1/T * len(G.graph['communities'][C])/len(G.nodes()) * sum(y[(t, u)] 
                                                                                    for u in G.nodes()) 
                        for t in range(T)) >= -0.1)
    
    m.setObjective(1/T * sum(sum(y[(t, v)] for t in range(T))
                            for v in G.nodes()), gp.GRB.MAXIMIZE)
    m.optimize()
    
    return [v for v in G.nodes() if x[v].X > 0]
    
