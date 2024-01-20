# modified from https://github.com/arnaucube/math/blob/master/bls12-381.sage
# and https://github.com/darkrenaissance/darkfi/blob/master/script/research/ec/pairing-modified-tate.sage
# and https://github.com/rdubois-crypto/FreshCryptoLib/blob/master/sage/FCL_pairings/arithmetic/curves/bls12_381.sage

from sage.all import *

# Generator of order r in E1 / F1
G1x = 0x17F1D3A73197D7942695638C4FA9AC0FC3688C4F9774B905A14E3A3F171BAC586C55E83FF97A1AEFFB3AF00ADB22C6BB
G1y = 0x8B3F481E3AAA0F1A09E30ED741D8AE4FCF5E095D5D00AF600DB18CB2C04B3EDD03CC744A2888AE40CAA232946C5E7E1

# Generator of order r in E2 / F2
G2x0 = 0x24AA2B2F08F0A91260805272DC51051C6E47AD4FA403B02B4510B647AE3D1770BAC0326A805BBEFD48056C8C121BDB8
G2x1 = 0x13E02B6052719F607DACD3A088274F65596BD0D09920B61AB5DA61BBDC7F5049334CF11213945D57E5AC7D055D042B7E
G2y0 = 0xCE5D527727D6E118CC9CDC6DA2E351AADFD9BAA8CBDD3A76D429A695160D12C923AC9CC3BACA289E193548608B82801
G2y1 = 0x606C4A02EA734CC32ACD2B02BC28B99CB3E287E85A763AF267492AB572E99AB3F370D275CEC1DA1AAA9075FF05F79BE

class BLS12381:
    def __init__(self) -> None:
        # self.x = -0xd201000000010000
        self.p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
        self.r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
        self.h1 = 0x396c8c005555e1568c00aaab0000aaab
        self.h2 = 0x5d543a95414e7f1091d50792876a202cd91de4547085abaa68a205b2e5a7ddfa628f1cb4d9e82ef21537e293a6691ae1616ec6e786f0c70cf1c38e31c7238e5

        self.F1 = GF(self.p)
        F2.<u> = GF(self.p^2, x, x^2+1)
        self.u = u
        self.F2 = F2
        F12.<w> = GF(self.p^12, x, x^12-2*x^6+2)
        # R12.<y> = PolynomialRing(F2)
        # F12.<w> = F2.extension(y^6-(u+1))
        self.w = w
        self.F12 = F12

class Pairing:
    def __init__(self, curve: BLS12381) -> None:
        self.curve = curve

        self.E1 = EllipticCurve(curve.F1, [0, 4])
        self.E2 = EllipticCurve(curve.F2, [0, 4*(1+curve.u)])
        self.E12 = EllipticCurve(curve.F12, [0, 4])
        self.E1.random_pointt = lambda: find_random_point(self.E1, self.curve.F1, 0, 4)
        self.E2.random_pointt = lambda: find_random_point(self.E2, self.curve.F2, 0, 4*(self.curve.u + 1))
        self.G1 = self.E1(G1x, G1y)
        self.G2 = self.E2(G2x0+curve.u*G2x1, G2y0+curve.u*G2y1)

    def lift_E1_to_E12(self, P):
        """
        Lift point on `E/F_q` to `E/F_{q^12}` using the natural lift
        """
        assert P.curve() == self.E1, "Attempting to lift a point from the wrong curve."
        return self.E12(P)

    def lift_E2_to_E12(self, P):
        """
        Lift point on `E/F_{q^2}` to `E/F_{q^12}` through the sextic twist
        """
        assert P.curve() == self.E2, "Attempting to lift a point from the wrong curve."
        # print(c.polynomial().coefficients() for c in (self.curve.h2*P).xy())
        xs, ys = [c.polynomial().coefficients() for c in (self.curve.h2*P).xy()]
        nx = self.curve.F12(xs[0] - xs[1] + self.curve.w ^ 6*xs[1])
        ny = self.curve.F12(ys[0] - ys[1] + self.curve.w ^ 6*ys[1])
        return self.E12(nx / (self.curve.w ^ 2), ny / (self.curve.w ^ 3))

    def pair(self, A, B):
        """
        Usage:

        ```
        bls12381 = BLS12381()
        e = Pairing(bls12381)
        r1 = e.E1.random_point()*bls12381.h1
        r2 = e.E2.random_point()*bls12381.h2
        r12 = e.pair(r1, r2)
        ```
        """
        A = self.lift_E1_to_E12(A)
        B = self.lift_E2_to_E12(B)
        return A.ate_pairing(B, self.curve.r, 12, self.E12.trace_of_frobenius())


def find_random_point(E, F, A, B):
    while True:
        x = F.random_element()
        y = sqrt(x ^ 3 + A * x + B)
        return E(x, y)
