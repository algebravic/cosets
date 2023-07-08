"""
Write a networkx undirected graph to DIMACS format
"""
from typing import Iterable, Any
from pathlib import Path
import networkx as nx

def _gen_dimacs(gph: nx.Graph) -> Iterable[str]:
    """
    Generate the lines of a DIMACS graph.
    """
    # Make a dictionary for the nodes
    node_map = {elt: ind for ind, elt in enumerate(gph.nodes, start=1)}
    yield from (f'c {elt}: {ind}' for elt, ind in node_map.items())
    yield f'p edge {len(gph.nodes)} {len(gph.edges)}'
    yield from (f'e {node_map[node1]} {node_map[node2]}'
                for node1, node2 in gph.edges)

def write_dimacs(gph: nx.Graph, name: str):
    """
    Write a DIMACS graph
    """
    with open(name, 'w', encoding='utf8') as fil:
        fil.write('\n'.join(_gen_dimacs(gph)))
        fil.write('\n')

def _gen_metis(gph: nx.Graph) -> Iterable[str]:
    """
    Generate the lines for a METIS graph.
    """
    ngph = nx.convert_node_labels_to_integers(gph, first_label=1)
    yield f'{len(ngph.nodes)} {len(ngph.edges)}'
    yield from (' '.join(map(str,sorted(ngph.neighbors(_)))) for _ in sorted(ngph.nodes))

def write_metis(gph: nx.Graph, name: str):
    """
    Write a METIS graph.
    """
    nfile = Path(name)
    with open(nfile.parent / (nfile.stem + '.metis'), 'w', encoding='utf8') as fil:
        fil.write('\n'.join(_gen_metis(gph)))
        fil.write('\n')
        
def normalize(elt: Any) -> str:
    """ Normalize list and tuple """
    return ('|'.join(map(str, elt))
            if isinstance(elt, (tuple, list))
            else str(elt))

def generate_csv(gph: nx.Graph) -> Iterable[str]:
    """ Generate the lines for a undirected graph """
    # Find isolated nodes
    # for node in gph.nodes:
    #     if gph.degree(node) == 0:
    #         yield normalize(node)
    for edge in gph.edges:
        yield normalize(edge[0]) + ',' + normalize(edge[1])

def write_csv(gph: nx.Graph, name: str):
    """
    Write a csv file for a graph.
    """
    with open(name, 'w', encoding='utf8') as fil:
        fil.write('\n'.join(generate_csv(gph)))
        fil.write('\n')

def write_independent(num: int, direct: str):
    """
    Write a simple LAD file for an independent graph """
    with open(direct + f'/i{num}.lad', 'w', encoding='utf8') as fil:
        fil.write(f'{num}\n')
        fil.write('\n'.join(num * ['0']))
        fil.write('\n')
