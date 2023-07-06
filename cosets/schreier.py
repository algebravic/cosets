"""Use Schreier-Sims to do full symmetry breaking.

After requiring 0 to be in the independent set,
the stabilizer of 0 is is isomorphic to S_2 x S_{n-1},
where the action is on the individual coordinates.
The generators that we'll use are those for S_n
as a Coxeter group (with 0-origin)
{(0,1),(2,3),(3,4), ..., (n-1,n)}
We'll then compute their induced action on the nodes of the graph,
and then give them to the Schreier-Sims algorithm.

Strategy:
1) We know we can include 0 in the independent set.
2) Inductively, if we've include, x[0], x[1], ...
in the independent set, we let G the stabilizer of these.
We then calculate the orbits of nodes in the graph which
are *not neighbors* of any of the x's, and then, inductively,
choose an orbit representative.  We thus form a tree.
In practice we can cut off the tree if it gets too big.

More details:

We are given a finite connected undirected graph G,
and a group of automorphisms, H, of G.
We want to find a maximum indpendent set.

Inductively, we have a sequence x[1], ..., x[m] of distinct nodes of G
which are an independent set.  Let H[m] denote the subgroup of H which fixes x[1], ..., x[m] and the graph G[m] := G \ {x[1], ..., x[m]}.
Let the orbits of H[m] acting on G[m] be O[1], ..., O[r], and
a[1] in O[1], ..., a[r] in O[r] be a choice of elements.  Then
we can choose x[m+1] to be one of a[1], ..., a[r].

Ab inito:

Proposition:
Let G be a connected undirected finite graph, and H its automorphism group.  Let the orbits of H acting on V(G) be O_1, ..., O_m, and
v_i in O_i be a representative.  Then a maximum independent subset of G
is H-equivalent to one which contains one of the elements v_i.

Proof: Let S denote a maximum independent set of G, and x in S.  By definition of orbit, there is an i such that x in O_i.  Let h in H be such that h x = v_i.  Then the set hS must contain v_i.

Comment: Let v in G, and let H be the stabilizer of v in Aut(G).  Then
the neighbors of v are a union of H orbits.

Proof: Since H is a subgroup of Aut(G) it permutes the edges of G and
nonedges of G.  If w is a neighbor of v and u a non-neighbor by H
cannot send the edge connecting (v,w) to the nonedge (v,u).

Corollary: Let G' denote the subgraph of G obtained by removing v and
all of its neighbors.  Then H stabilizes G'

Proof: The action of H on the nodes G forms orbits.  We have already
seen that two of the orbits are {v} and {w: w is a neighbor of v}.
Thus the nodes of G' are a union of orbits of H.

Lemma: Let G be a graph whose nodes are labeled by all the
n-tuples of 0/1, and S be a set of n-tuples.  The edges of G
are all of the form (x, x ^ s) for all x and and all s in S
where ^ denotes elementwise XOR.  Let H denote the subgroup
of S_n which leaves S invariant.  Then H is a subgroup of
Aut(G).

Proof: Let h in H, and (x, x ^ s) be an edge of G.
Then (h(x), h (x ^ s)) = (h(x), h(x) ^ h(s)) which
is an edge of G since h(s) in S.

Question: Can there be any other automorphisms?
If g is an automorphism we must have
(g(x), g(x ^ s)) is an edge for all x, s.
That is g(x) ^ g(x ^ s)  in S. I guess that this could happen.

Plan: Recursively find the orbits of the graph.
Choose a representative of each orbit.
For each such recursively find its stabilizer, and then the
reduced graph obtained by deleting it and its neighbors, and
then the reduced automorphism group of the reduced graph.
[Probably an easily described subgroup].
Choose an element (probably the least).

New Plan: Our goal is to produce a rooted tree, where
each node, except the root, is labeled by a distinct node in the
graph.  The tree has the following property:
Any maximum size independent set is equivalent under the
automorphism group of the graph to one with the following property:
For every path from the root the leaf, conditional on all the
internal nodes in the path the leaf is included.

We will construct the tree as follows:

1) Calculate the automorphism group of the graph.
2) For orbit of the nodes under the automorphism group
choose a representative.
3) Recursively, for each representative, create a new graph
from the old by removing that node and all its neighbors,
and go back to (1).

We'll need a practical termination criterion.
Maybe a cutoff parameter.

Our goal is to generate a tree.  We assume that the input
graph is vertex transitive (for now).  Each node of the tree
will contain the following information
1) A graph
2) A node
3) An automorphism group, the stabilizer of the node.
4) The cardinality of the orbit of the node under the automorphism group

A node will have the following children:
a) if the cardinality of the orbit is 1, it is a leaf.
b) If greater than 1, the following children are genearted:
first generate a new graph, by removing all of the neighbors
of the node, and the node, itself from the parent graph.
For each orbit of the action of the automorphism group on the
new graph, choose a representative (such as the minimum element),
and generate a new automorphism group which is the stabilizer
of the parent automorphism group acting on the representative node.
n
Note: We would like to have a new group which is not the stabilizer
of the node, but the stabilizer of the subset of all nodes
from this tree node to the root.  However, this appears to be
too onerous a calculation.

Expanding the tree completely is generally not feasible, so
we will have truncation rules.

For any truncation of the tree, the following is true:
If A is a maximum cardinality indepdendent set, there is some
element g in the root automorphism group so that A^g contains
all of the nodes in some path from a leaf to the root.

"""
from typing import Iterable, List, Tuple, Set, Callable
from functools import partial
from itertools import product, chain
from collections import namedtuple
import networkx as nx
from pysat.formula import IDPool
from lazytree import LazyTree
from sympy.combinatorics import Permutation, PermutationGroup

POINT = Tuple[int, ...]
CLAUSE = List[int]
TreeNode = namedtuple('TreeNode',
                      ['node', 'number', 'graph', 'group'])

def subset_stabilizer(grp: PermutationGroup,
                      points: Set[int]) -> PermutationGroup:
    """
    Find the stabilizer of the subset points by the group grp.
    """
    base, gens = grp.schreier_sims_incremental()
    return grp.subgroup_search(lambda elt:
                               {elt(_) for _ in points} == points,
                               base = base,
                               strong_gens = gens)

def choices(grp: PermutationGroup,
            support: Set[int]) -> Iterable[Tuple[int, int]]:
    """
    Reperesentatives of the orbits.
    """
    orbs = grp.orbits()
    for orb in orbs:
        if support.issuperset(orb):
            yield min(orb), len(orb)

def dgraph(graph: nx.Graph, node: int) -> nx.Graph:
    nbrs = set(list(graph.neighbors(node)))
    ngraph = graph.copy()
    for nbr in nbrs:
        ngraph.remove_node(nbr)
    ngraph.remove_node(node)
    return ngraph

def children(tnode: LazyTree) -> List[LazyTree]:
    """ Form the children. """

    clist = choices(tnode.group, set(list(tnode.graph.nodes)))
    return [TreeNode(graph = dgraph(tnode.graph, cnode),
                     group = tnode.group.stabilizer(cnode),
                     node = cnode,
                     number = cnum)
            for cnode, cnum in clist]

def make_tree(gph: nx.Graph,
              grp: PermutationGroup) -> LazyTree:
    """
    Make the stabilizer tree
    gph: the root graph.
    grp: a group of automorphism of gph for which gph is
         vertex transitive.
    """
    node = min(gph.nodes)
    return LazyTree(root = TreeNode(graph = dgraph(gph, node),
                                    group = grp.stabilizer(node),
                                    node = node,
                                    number = len(gph.nodes)),
                    child_map = children,
                    view = lambda _: _.node)

def tree_clauses(pool: IDPool,
                 test: Callable[[LazyTree], bool],
                 depth: int,
                 tree: LazyTree) -> Iterable[CLAUSE]:
    """
    Produce clauses for symmetry breaking from the tree, using
    a breadth first search
    """
    lit = pool.id(('x', tree.root.node))
    if depth > 0 and test(tree):
        below = list(tree.children) # don't recalculate
        yield ([-lit]
               + [pool.id(('x', child.root.node))
                  for child in below])
        for cls in chain(*(tree_clauses(pool, test, depth - 1, child)
                           for child in below)):
            yield [-lit] + cls
    
def transposition(num: int, inds: Tuple[int, int]) -> List[int]:
    """
    Transposition.
    """
    ind, jind = inds
    val = list(range(num))
    val[ind], val[jind] = jind, ind
    return val
