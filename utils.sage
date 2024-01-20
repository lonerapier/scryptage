from sage.all import *


def random_curve_values(E, n):
    return [E.random_element() for i in range(n)]


def random_field_values(F, n):
    return [F.random_element() for i in range(n)]


def inner_product_vec(a, b):
    assert len(a) == len(b)

    res = 0
    for i in range(len(a)):
        res += a[i] * b[i]
    return res


def inner_product_point(a, b):
    assert len(a) == len(b)
    c = 0
    for i in range(len(a)):
        c = c + (a[i]) * b[i]
    return c


def scalar_mul_field(f, v, n):
    res = [None] * n
    for i in range(n):
        res[i] = f[i] * v
    return res


def scalar_mul_point(g, v, n):
    r = [None] * len(g)
    for i in range(len(g)):
        r[i] = g[i] * (v)
    return r


def add_vectors(a, b):
    assert len(a) == len(b)
    return [x + y for x, y in zip(a, b)]


def powers_of(x, n):
    return [x ^ i for i in range(n)]
