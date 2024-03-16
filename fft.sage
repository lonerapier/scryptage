# https://www.csd.uwo.ca/~mmorenom/CS874/Lectures/Newton2Hensel.html/node9.html

from sage.all import *
import unittest


def generator(F):
    q = F.order()
    factors = list(factor(q - 1))
    for i in range(2, q):
        g = F(i)

        # Check if g is a generator by checking if g ^ ((q-1)/f) != 1 for all factors f

        flag = True
        for f in factors:
            if g ^ ((q - 1) / f[0]) == 1:
                flag = False

        if flag:
            return g


def primitive_rou(F, n):
    g = generator(F)
    # print(g)
    return g ^ ((F.order() - 1) / n)


def vandermonde_matrix(F, n, omega):
    return matrix(F, n, n, lambda i, j: omega ^ (i * j))


def vandermonde_det(M):
    prod = 1
    for i in range(M.nrows()):
        for j in range(i + 1, M.nrows()):
            # print(M[i, 1], M[j, 1])
            prod *= M[i, 1] - M[j, 1]

    # print(f"prod: {prod}")
    return prod


def fft1(F, n, omega, a):
    """
    calculates the fft of a polynomial using the vandermonde matrix
    inputs:
    - F: finite field
    - n: degree of the polynomial
    - omega: primitive n-th root of unity
    - a: polynomial
    """

    M = vandermonde_matrix(F, n, omega)

    return M * vector(a)


def ifft1(F, n, omega, a):
    """
    calculates the ifft of a polynomial using the vandermonde matrix
    inputs:
    - F: finite field
    - n: degree of the polynomial
    - omega: primitive n-th root of unity
    - a: evaluation of the polynomial at the n-th roots of unity
    """

    M = vandermonde_matrix(F, n, omega)

    return M.inverse() * vector(a)


def fft2(F, n, omega, a):
    """
    calculates the fft of a polynomial using the cooley-tukey algorithm

    inputs:
    - F: finite field
    - n: degree of the polynomial
    - omega: primitive n-th root of unity
    - a: polynomial
    """

    if n == 1:
        return a

    n2 = n // 2
    omega2 = omega ^ 2

    a_even = [a[i] for i in range(0, n, 2)]
    a_odd = [a[i] for i in range(1, n, 2)]

    y_even = fft2(F, n2, omega2, a_even)
    y_odd = fft2(F, n2, omega2, a_odd)

    y = [F(0) for _ in range(n)]
    for k in range(n2):
        y[k] = y_even[k] + omega ^ k * y_odd[k]
        y[k + n2] = y_even[k] - omega ^ k * y_odd[k]

    return y


def sample_poly(F):
    a = [
        F(1),
        F(2),
        F(3),
        F(4),
        F(5),
        F(6),
        F(7),
        F(8),
        F(9),
        F(10),
        F(11),
        F(12),
        F(13),
        F(14),
        F(15),
        F(16),
    ]
    return a


class TestFFT(unittest.TestCase):
    def test_generator(self):
        p = 17
        F = GF(p)
        g = generator(F)
        self.assertEqual(g.multiplicative_order(), 16)

        F = GF(2 ^ 64 - 2 ^ 32 + 1)
        g = generator(F)
        self.assertEqual(g.multiplicative_order(), 2 ^ 64 - 2 ^ 32)

    def test_primitive_rou(self):
        p = 17
        F = GF(p)
        g = primitive_rou(F, 16)
        self.assertEqual(g.multiplicative_order(), 16)

        F = GF(2 ^ 64 - 2 ^ 32 + 1)
        g = primitive_rou(F, 2 ^ 32)
        self.assertEqual(g.multiplicative_order(), 2 ^ 32)

    def test_vandermonde_matrix(self):
        p = 17
        F = GF(p)
        g = primitive_rou(F, 16)
        M = vandermonde_matrix(F, 16, g)
        self.assertEqual(M.det(), vandermonde_det(M))

        F = GF(2 ^ 64 - 2 ^ 32 + 1)
        g = primitive_rou(F, 2 ^ 32)
        # M = vandermonde_matrix(F, 2 ^ 32, g)
        # self.assertEqual(M.det(), vandermonde_det(M))

    def test_fft1(self):
        p = 17
        F = GF(p)
        g = primitive_rou(F, p - 1)
        a = sample_poly(F)
        R = PolynomialRing(F, "x")
        a_poly = R(a)
        a_eval = fft1(F, p - 1, g, a)

        print(a_eval)
        for i in range(p - 1):
            self.assertEqual(a_eval[i], a_poly(g ^ i))

    def test_ifft1(self):
        p = 17
        F = GF(p)
        g = primitive_rou(F, p - 1)
        a = sample_poly(F)
        R = PolynomialRing(F, "x")
        a_poly = R(a)
        a_eval = fft1(F, p - 1, g, a)
        a_recover = ifft1(F, p - 1, g, a_eval)

        print(a_recover)
        for i in range(p - 1):
            self.assertEqual(a_recover[i], a[i])

    def test_fft2(self):
        p = 17
        F = GF(p)
        a = sample_poly(F)
        g = primitive_rou(F, p - 1)
        R = PolynomialRing(F, "x")
        a_poly = R(a)
        a_eval = fft2(F, p - 1, g, a)

        for i in range(p - 1):
            self.assertEqual(a_eval[i], a_poly(g ^ i))


if __name__ == "__main__":
    unittest.main()
