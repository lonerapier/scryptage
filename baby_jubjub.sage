# modified from https://github.com/bellesmarta/baby_jubjub/blob/master/findCurve.sage
# and https://github.com/iden3/iden3-docs/blob/master/source/iden3_repos/research/publications/zkproof-standards-workshop-2/baby-jubjub/Baby-Jubjub.pdf

P = 21888242871839275222246405745257275088548364400416034343698204186575808495617
Fr = GF(P)
h = 8  # cofactor

A = 168698
EC = EllipticCurve(Fr, [0, A, 0, 1, 0])

DTed = 168696
ATed = 168700


N = EC.order()
assert N % h == 0
print("order", N, "\nfactor:", factor(N))


def NewPoint(curve):
    return curve(0, 1)


def FindGenPoint(F, EC, A, N):
    i = 1
    while True:
        u = F(i)
        v2 = u ^ 3 + A * u ^ 2 + u
        if not v2.is_square():
            i = i + 1
            continue
        v = v2.sqrt()
        point = EC(u, v)
        i = i + 1
        if point.order() == N:
            return point


xGen, yGen, _ = FindGenPoint(Fr, EC, A, N)


def MontToTEdw(u, v, r):
    x = Mod(u / v, r)
    y = Mod((u - 1) / (u + 1), r)
    return (x, y)


def TedwToMont(x, y, r):
    u = Mod((1 + y) / (1 - y), r)
    v = Mod((1 + y) / ((1 - y) * x), r)
    return (u, v)


def IsOnTEdw(x, y, r, a, d):
    return (
        Mod(
            Mod(a, r) * Mod(x ^ 2, r)
            + Mod(y ^ 2, r)
            - 1
            - Mod(d, r) * Mod(x ^ 2, r) * Mod(y ^ 2, r),
            r,
        )
        == 0
    )


xGenTed, yGenTed = MontToTEdw(xGen, yGen, P)

assert IsOnTEdw(xGenTed, yGenTed, P, ATed, DTed) == True


# EC2 = EC.quadratic_twist()
# print(EC2)
