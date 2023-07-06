"""
Calculate the Nemhauser-Trotter relaxation for independent set.
"""
from typing import Dict, Hashable
import networkx as nx
from networkx import bipartite

def nt_relax(gph: nx.Graph) -> Dict[Hashable, int]:
    """
    Form the bipartite graph, whose right and left
    vertices correspond to the nodes in the graph
    and are connected with an edge if they are connected
    in the original graph.
    """

    bgph = nx.Graph()
    for node in gph.nodes:
        bgph.add_node(('l', node), bipartite=0)
        bgph.add_node(('r', node), bipartite=1)
    for node1, node2 in gph.edges:
        bgph.add_edge(('l', node1), ('r', node2))
        bgph.add_edge(('r', node1), ('l', node2))

    matching = bipartite.maximum_matching(bgph)
    cover = bipartite.to_vertex_cover(bgph, matching)
    out = {}
    for node in gph.nodes:

        tst = set([('l', node), ('r', node)])
        out[node] = len(tst.intersection(cover))
    return out
