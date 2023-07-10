"""
Use the mip package to model independent set
"""
from typing import Iterable, Hashable, Tuple, Dict
import networkx as nx
from mip import Model, INTEGER, CONTINUOUS, BINARY, xsum, MAXIMIZE

def mip_model(gph: nx.Graph) -> Tuple[Model, Dict[Hashable, int]]:
    """
    Use the standard mip_model.
    """
    ngph = nx.convert_node_labels_to_integers(gph,
                                              ordering='sorted')
    dct = {ind: val for ind, val in enumerate(sorted(gph.nodes))}
    mvars = {}
    model = Model(sense=MAXIMIZE)
    mvars = [model.add_var(name = f'x{node}', var_type=BINARY)
             for node in ngph.nodes]
        
    for node1, node2 in ngph.edges:
        model += mvars[node1] + mvars[node2] <= 1
    model.objective = xsum(mvars)

    return model, dct

def cnf_model(cnf: WCNF) -> Tuple[Model, Dict[int, int]]:
    """
    Render a weight CNF as a MIP
    """
    top = cnf.nv
