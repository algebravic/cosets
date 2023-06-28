"""
The distance graph of Z^n/2D_n.
"""
from typing import Tuple, Iterable, Any
from itertools import chain, product
from sympy import binomial

VEC = Tuple[int,...]

def small_weight(num: int, wgt: int) -> Iterable[VEC]:
    """
    All binary vectors of weight wgt of length n.
    """
    if wgt == 0:
        yield num * (0,)
    elif wgt == num:
        yield num * (1,)
    elif 0 < wgt < num:
        for elt in range(0,2):
            for rest in small_weight(num - 1, wgt - elt):
                yield (elt,) + rest

def dn_neighbors(num: int) -> Iterable[VEC]:
    """
    Neighbors of 0 in the Dn graph    
    """

    yield from (delta + (num - 1) * (0,)
                for delta in [(0,1), (1,0)])
    yield from (pref + delta
                for pref in product(range(2), repeat=2)
                for delta in chain(small_weight(num - 1, 1),
                                   small_weight(num - 1, 2)))
    yield from (ndelta + delta
                for ndelta in [(0,0), (1,1)]
                for delta in small_weight(num - 1, 3))

def dn_counting(num: int) -> int:
    """
    Count the number of edges per Veit Elser.
    """
    cnt = 2
    cnt += 4 * binomial(num-1,1)
    cnt += 4 * binomial(num-1,2)
    cnt += 2 * binomial(num-1,3)
    return (2 ** num) * cnt
