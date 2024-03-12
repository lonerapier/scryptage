from sage.rings.polynomial.polydict import PolyDict


def multiplicative_identity(num):
    return PolyDict({tuple([0] * num): 1})


def additive_identity(num):
    return PolyDict({tuple([0] * num): 0})
