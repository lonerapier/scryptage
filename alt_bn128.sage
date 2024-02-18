# modified from https://github.com/scipr-lab/libff/blob/develop/libff/algebra/curves/alt_bn128/alt_bn128.sage
# and https://github.com/ethereum/py_pairing/blob/master/py_ecc/bn128/bn128_field_elements.py
from sage.all import *

G1x = 1
G1y = 2

G2x0 = 10857046999023057135944570762232829481370756359578518086990519993285655852781
G2x1 = 11559732032986387107991004021392285783925812861821192530917403151452391805634
G2y0 = 8495653923123431417604973247489272438418190587263600148770280649306958101930
G2y1 = 4082367875863433681332203403145435568316851327593401208105741076214120093531


def g2_h(x):
    return 36 * x ^ 4 + 36 * x ^ 3 + 30 * x ^ 2 + 6 * x + 1


class BN254:
    """
    ## Usage:
    ```
    bn254 = BN254()
    print(bn254.p, bn254.x)
    ```
    """

    def __init__(self) -> None:
        # curve parameter
        self.x = 4965661367192848881
        # curve prime
        self.p = 21888242871839275222246405745257275088696311157297823662689037894645226208583
        # trace of frobenius
        self.t = 29793968203157093285
        # curve order
        self.r = 21888242871839275222246405745257275088548364400416034343698204186575808495617

        self.Fp = GF(self.p)
        self.Fr = GF(self.r)
        # F2.<u> = GF(self.p^2, x, x^2+1)
        # self.u = u
        # self.F2 = F2
        # F12.<w> = GF(self.p^12, x, x^12-18*x^6+82)
        # self.w = w
        # self.F12 = F12

        self.E = EllipticCurve(self.Fp, [0, 3])
        # self.E2 = EllipticCurve(self.F2, [0, 3/(9+self.u)])
        # self.E12 = EllipticCurve(self.F12, [0, 3])

        N = self.E.order()

        self.h1 = 1
        self.h2 = g2_h(self.x)

        self.G1 = self.E(G1x, G1y)
        assert self.G1.curve() == self.E
        assert self.G1.order() == self.r

        # self.G2 = self.E2(G2x0+self.u*G2x1, G2y0+self.u*G2y1)

        # assert self.G2.curve() == self.E2

        self.check_p(self.x)
        # self.check_t(self.x)
        self.check_r(self.x)
        self.check_g1_cofactor()

    def check_p(self, x):
        p = 36 * (x ^ 4) + 36 * (x ^ 3) + 24 * (x ^ 2) + 6 * x + 1
        assert self.p == p

    def check_t(self, x):
        t = 6 * x ^ 2 + 1
        assert self.t == t

    def check_r(self, x):
        r = 36 * x ^ 4 + 36 * x ^ 3 + 18 * x ^ 2 + 6 * x + 1
        assert self.r == r

    def check_g1_cofactor(self):
        order = factor(self.r)
        assert order[0] == (self.r, 1)

    # def twist_E2_E12(self, p):
    #     assert p.curve() == self.E2


bn254 = BN254()
print(bn254.p, bn254.x)
