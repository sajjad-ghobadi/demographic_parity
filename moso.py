
import numpy as np

from tsang.algorithms import make_normalized, maxmin_algo
from tsang.icm import (make_multilinear_gradient_group,
                       make_multilinear_objective_samples_group)


def tsang_maximin_md(G):
    """Call Tsang et al 2019's implementation on G using md as solver.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.

    Returns
    -------
    dict
        dictionary with keys nodes and values probabilities

    """
    return moso(G, solver='md')


def moso(G, solver='gurobi'):
    """Call Tsang et al 2019's implementation on G using solver as solver.

    Parameters
    ----------
    G : networkx.DiGraph
        The underlying considered instance.
    solver: str
        string encoding the solver to be used, either gurobi or md

    Returns
    -------
    dict
        dictionary with keys nodes and values probabilities

    """
    budget = G.graph['k']
    batch_size = 1000
    threshold = 5

    live_graphs = G.graph['leg_tsang']

    group_indicator = np.zeros((len(G.nodes()), len(G.graph['communities'])))
    community_i = 0
    for community in G.graph['communities'].keys():
        for node in G.graph['communities'][community]:
            group_indicator[node, community_i] = 1
        community_i += 1
    community_sizes = sum(group_indicator)

    val_oracle = make_multilinear_objective_samples_group(
        live_graphs, group_indicator,  sorted(list(G.nodes())), sorted(list(G.nodes())),
        np.ones(len(G)))
    grad_oracle = make_multilinear_gradient_group(
        live_graphs, group_indicator,  sorted(list(G.nodes())), sorted(list(G.nodes())),
        np.ones(len(G)))

    grad_oracle_normalized = make_normalized(grad_oracle, community_sizes)
    val_oracle_normalized = make_normalized(val_oracle, community_sizes)
    minmax_x = maxmin_algo(grad_oracle_normalized, val_oracle_normalized,
                           threshold, budget, group_indicator, 20, 10, 0.05,
                           solver, batch_size)

    minmax_x = minmax_x.mean(axis=0)

    return {i: minmax_x[i] for i in G.nodes()}
