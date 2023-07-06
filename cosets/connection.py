"""
The distance graph of Z^n/2D_n.
"""
from typing import Tuple, Iterable, Any, List
from itertools import chain, product
from sympy import binomial

def candidates(num: int) -> List[VEC]:
    """
    Candidates for the equivalence class after removing 0.
    Brute force
    """
    good = set()
    universe = set(product(range(2), repeat = num + 1))
    consider = universe.difference(list(dn_neighbors(num))
                                   + [(num + 1) * (0,)])
    # representative after permuting first 2 coords and last n-1
    for elt in consider:
        equiv = tuple(list(sorted(elt[:2]))
                      + list(sorted(elt[2:])))
        good.add(equiv)
    return list(sorted(good))
