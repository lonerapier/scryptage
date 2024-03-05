# modified from https://github.com/arnaucube/math/blob/master/ipa.sage
from sage.all import *

load("alt_bn128.sage")
load("utils.sage")


def construct_s_from_u(u, d):
    k = int(log(d, 2))
    s = [1] * d

    t = d
    for i in reversed(range(k)):
        t /= 2
        x = 0
        for j in range(d):
            if x < t:
                s[j] *= u[i] ^ -1
            else:
                s[j] *= u[i]
            x = x + 1
            if x >= t * 2:
                x = 0
    return s


def amortization_poly(F, u):
    R.<x> = PolynomialRing(F)
    for i in range(1, len(u)+1):
        u_k = u[i-1]
        u_k_inv = u[i-1] ^ -1
        R *= ((u_k_inv * (x ** (2 ^ (i - 1)))) + u_k)
    return R

class Halo:
    def __init__(self, F, E, N) -> None:
        self.F = F
        self.E = E
        self.N = N

    def setup(self):
        k = int(log(self.N, 2))
        self.G = random_curve_values(self.E, self.N)
        self.H = self.E.random_element()
        self.U = self.E.random_element()
        u = random_field_values(self.F, k)
        return u

    def inner_product_relation(self, a, b):
        r = self.F.random_element()
        P = (
            inner_product_point(a, self.G)
            + r * self.H
            + (inner_product_vec(a, b)) * self.U
        )
        return P, r

    def ipa(self, _a, _b, u):
        a = _a
        b = _b
        g = self.G

        k = int(log(self.N, 2))

        l = [None] * k
        r = [None] * k
        L = [None] * k
        R = [None] * k

        for i in reversed(range(0, k)):
            print("i: {}, n: {}, m: {}".format(i, len(a), len(a)/2))
            m = int(len(a) / 2)
            a_l = a[:m]
            a_r = a[m:]
            b_l = b[:m]
            b_r = b[m:]
            g_l = g[:m]
            g_r = g[m:]

            l[i] = self.F.random_element()
            r[i] = self.F.random_element()

            L[i] = (
                inner_product_point(a_l, g_r)
                + l[i] * self.H
                + (inner_product_vec(a_l, b_r)) * self.U
            )
            R[i] = (
                inner_product_point(a_r, g_l)
                + r[i] * self.H
                + (inner_product_vec(a_r, b_l)) * self.U
            )

            u_k = u[i]
            u_k_inv = u[i] ** (-1)

            a = add_vectors(
                scalar_mul_field(a_l, u_k, m), scalar_mul_field(a_r, u_k_inv, m)
            )
            b = add_vectors(
                scalar_mul_field(b_l, u_k_inv, m), scalar_mul_field(b_r, u_k, m)
            )
            g = add_vectors(
                scalar_mul_point(g_l, u_k_inv, m), scalar_mul_point(g_r, u_k, m)
            )

        assert len(a) == 1
        assert len(b) == 1
        assert len(g) == 1
        assert len(r) == k

        return a[0], b[0], g[0], l, r, L, R

    def verify(self, P, u, a, _b, powers_of_x, g, r, l_prime, r_prime, L, R):
        LHS = P
        s = construct_s_from_u(u, self.N)
        b = inner_product_vec(powers_of_x, s)
        G = inner_product_point(s, self.G)

        assert(b == _b)

        # amort_poly = amortization_poly(self.F, u)
        # b_amort = amort_poly
        # print("amort_poly: {}\nb_amort: {}".format(amort_poly, b_amort))

        r_dash = r
        RHS = P
        for i in range(len(L)):
            u_k = u[i]^2
            u_k_inv = u[i] ^ -2
            r_dash += l_prime[i]*u_k + r_prime[i]*u_k_inv
            RHS = RHS + (u_k * L[i]) + (u_k_inv * R[i])

        LHS = a*G + r_dash * self.H + a*b*self.U

        assert(LHS == RHS)

bn254 = BN254()
Fq = bn254.Fr

n = 8
a = [
    Fq(1),
    Fq(2),
    Fq(3),
    Fq(4),
    Fq(5),
    Fq(6),
    Fq(7),
    Fq(8),
]
x = Fq(3)
b = [x**i for i in range(n)]

halo_ipa = Halo(bn254.Fr, bn254.E, n)

u = halo_ipa.setup()

P, r = halo_ipa.inner_product_relation(a, b)

a_0, b_0, g_0, l_prime, r_prime, L, R = halo_ipa.ipa(a, b, u)

halo_ipa.verify(P, u, a_0, b_0, b, g_0, r, l_prime, r_prime, L, R)
