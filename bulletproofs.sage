from sage.all import *

load("utils.sage")


class Bulletproofs:
    def __init__(self, F, E, n) -> None:
        self.F = F
        self.E = E
        self.N = n

    def setup(self, g1):
        self.G = random_curve_values(self.E, self.N)
        self.H = random_curve_values(self.E, self.N)
        self.U = self.E.random_element()
        u = random_field_values(self.F, self.N)

        # temp for debugging purpose
        #
        # self.G = [g1, g1, g1, g1]
        # self.H = [g1, g1, g1, g1]
        # self.U = g1
        # u = [self.F(1), self.F(1), self.F(1), self.F(1)]

        return u

    def inner_product_relation(self, a, b):
        P = (
            inner_product_point(a, self.G)
            + inner_product_point(b, self.H)
            + (inner_product_vec(a, b)) * self.U
        )
        return P

    def ipa(self, _a, _b, u):
        a = _a
        b = _b
        g = self.G
        h = self.H

        k = log(self.N, 2)

        L = [None] * k
        R = [None] * k

        for i in reversed(range(0, k)):
            m = int(len(a) / 2)
            a_l = a[:m]
            a_r = a[m:]
            b_l = b[:m]
            b_r = b[m:]
            g_l = g[:m]
            g_r = g[m:]
            h_l = h[:m]
            h_r = h[m:]

            L[i] = (
                inner_product_point(a_l, g_r)
                + inner_product_point(b_r, h_l)
                + (inner_product_vec(a_l, b_r)) * self.U
            )
            R[i] = (
                inner_product_point(a_r, g_l)
                + inner_product_point(b_l, h_r)
                + (inner_product_vec(a_r, b_l)) * self.U
            )

            # print("i: {}, L[i]: {}, R[i]: {}".format(i, L[i], R[i]))

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
            h = add_vectors(
                scalar_mul_point(h_l, u_k, m), scalar_mul_point(h_r, u_k_inv, m)
            )

        assert len(a) == 1
        assert len(b) == 1
        assert len(g) == 1
        assert len(h) == 1

        return a[0], b[0], g[0], h[0], L, R

    def verify(self, P, u, a, b, g, h, L, R):
        LHS = P
        for i in range(len(L)):
            u_k = u[i]
            u_k_inv = u[i] ^ (-1)
            LHS = LHS + L[i] * (u_k**2) + R[i] * (u_k_inv**2)
        RHS = a * g + b * h + (a * b) * self.U

        print("LHS: {}\nRHS: {}".format(LHS, RHS))
        assert LHS == RHS


load("alt_bn128.sage")

bn254 = BN254()

n = 8
a = [
    bn254.Fr(1),
    bn254.Fr(2),
    bn254.Fr(3),
    bn254.Fr(4),
    bn254.Fr(5),
    bn254.Fr(6),
    bn254.Fr(7),
    bn254.Fr(8),
]
x = bn254.Fr(3)
b = [x**i for i in range(n)]

bulletproofs_ipa = Bulletproofs(bn254.Fr, bn254.E, n)

u = bulletproofs_ipa.setup(bn254.G1)

P = bulletproofs_ipa.inner_product_relation(a, b)

a_0, b_0, g_0, h_0, L, R = bulletproofs_ipa.ipa(a, b, u)

bulletproofs_ipa.verify(P, u, a_0, b_0, g_0, h_0, L, R)


class RangeProofs:
    def __init__(self) -> None:
        pass
