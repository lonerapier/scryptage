# https://github.com/darkrenaissance/darkfi/blob/master/script/research/zk/nova-fold.sage
# https://github.com/darkrenaissance/darkfi/blob/master/script/research/zk/nova-ivc.sage

from sage.all import *
import unittest

load("pasta.sage")


temp = F1
pallas = pallas

F1 = F2  # scalar field of E1
F2 = temp  # scalar field of E2
vesta = vesta

assert pallas.order() == F1
assert vesta.order() == F2

pallas_gen = pallas.gens()[0]
vesta_gen = vesta.gens()[0]


def hadamard_prod(F, a, b):
    assert len(a) == len(b)
    res = [a[i] * b[i] for i in range(len(a))]
    return vector(F, res)


def pedersen_commit(R, G, F, A):
    """
    pedersen vector committment

    inputs:
    - `R`: random points
    - `G`: generator point
    - `F`: field
    - `A`: vector to be committed
    """
    res = G * F.random_element()
    for i in range(len(A)):
        res += A[i] * R[i]
    return res


def pedersen_commit_table(table, R, G, F, A):
    if table[A] != "":
        return table[A]

    commit = pedersen_commit(R, F, G, A)
    table[A] = commit
    return commit


def evaluate(x, y, s):
    xy = x * y
    sxy = s * xy
    w1 = (1 - s) * (x + y)
    out = sxy + w1
    return out


# arithmetisation for:
# sxy + (1 - s)(x + y) = z
# s(1 - s) = 0
# W = [x, y, s, xy, xys, w1]
# X = [out]
# Z = (x, y, s, xy, xys, w1, out, 1)
def witness_instance_pair(F, x, y, s):
    x = F(x)
    y = F(y)
    s = F(s)

    xy = x * y
    sxy = s * xy
    w1 = (F(1) - s) * (x + y)
    out = sxy + w1

    W = [x, y, s, xy, sxy, w1]
    X = [out, F(1)]
    return W, X


# Circuit
A = Matrix(
    [
        [1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, -1, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 1, 1, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0],
    ],
)
B = Matrix(
    [
        [0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0],
        [1, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1],
        [0, 0, -1, 0, 0, 0, 0, 1],
    ],
)
C = Matrix(
    [
        [0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
    ],
)


class CommittedRelaxedR1CS:
    def __init__(self, comm_E, μ1, comm_W, x1):
        self.comm_E = comm_E
        self.μ = μ1
        self.comm_W = comm_W
        self.x = x1


class RelaxedR1CSWitness:
    def __init__(self, E, re, W, rw) -> None:
        self.E = E
        self.re = re
        self.W = W
        self.rw = rw


# section 4.2 of [Nova](https://eprint.iacr.org/2021/370)
class NIFS:
    def __init__(self, F, E, A, B, C) -> None:
        self.F = F
        self.E = E
        self.A = A
        self.B = B
        self.C = C
        self.R = [E.random_element() for _ in range(10)]

    def prover(self, Z1, Z2, r):
        """
        generate nova folding scheme proof for committed relaxed R1CS instance

        inputs:
        - `F`: field
        - `E`: curve for commit
        - `R`: random curve points for pedersen commit
        - `A`: R1CS variable A
        - `B`: R1CS variable B
        - `C`: R1CS output C
        - `Z1`: instance witness pair for computation 1
        - `Z2`: instance witness pair for computation 2

        outputs:
        - `E`: commitment to error term E
        - `T`: commitment to term T
        - `comm_w1`: commitment to W1
        - `comm_w2`: commitment to W2
        """
        μ1, μ2 = 1, 1
        e1, e2 = vector([0] * 5), vector([0] * 5)

        comm_E1 = pedersen_commit(self.R, pallas_gen, self.F, e1)
        comm_E2 = pedersen_commit(self.R, pallas_gen, self.F, e2)
        comm_w1 = pedersen_commit(self.R, pallas_gen, self.F, W1)
        comm_w2 = pedersen_commit(self.R, pallas_gen, self.F, W2)

        T = (
            hadamard_prod(self.F, A * Z1, B * Z2)
            + hadamard_prod(self.F, A * Z2, B * Z1)
            - μ1 * C * Z2
            - μ2 * C * Z1
        )
        comm_T = pedersen_commit(self.R, pallas_gen, self.F, T)

        E = e1 + r * T + (r ^ 2) * e2
        W = vector(W1) + r * vector(W2)

        r1 = CommittedRelaxedR1CS(comm_E1, μ1, comm_w1, X1)
        r2 = CommittedRelaxedR1CS(comm_E2, μ2, comm_w2, X2)
        return E, W, comm_T, r1, r2

    def verifier(self, comm_T, r, r1: CommittedRelaxedR1CS, r2: CommittedRelaxedR1CS):
        comm_W = r1.comm_W + r * r2.comm_W
        comm_E = r1.comm_E + r * comm_T + (r ^ 2) * r2.comm_E
        X = vector(r1.x) + r * vector(r2.x)
        μ = r1.μ + r * r2.μ

        return comm_E, μ, comm_W, X


def Func1(F, x):
    return x**5


def Func2(F, x):
    return x**3


hash_table1 = {}
hash_table2 = {}


def hash1(F, x):
    if hash_table1[x] != "":
        return hash_table1[x]

    val = F.random_element()
    while val > 2 ^ 250 - 1:
        val = F.random_element()

    hash_table1[x] = val
    return val


def hash2(F, x):
    if hash_table2[x] != "":
        return hash_table2[x]

    val = F.random_element()
    while val > 2 ^ 250 - 1:
        val = F.random_element()

    hash_table2[x] = val
    return val


commit_table1 = {}
commit_table2 = {}


def commit1(F, R, G, x):
    return pedersen_commit_table(commit_table1, R, G, F, x)


def commit2(F, R, G, x):
    return pedersen_commit_table(commit_table2, R, G, F, x)


class IVCProof:
    def __init__(self, u2: CommittedRelaxedR1CS, w2, U1, W1, U2, W2):
        self.u2 = u2
        self.w2 = w2
        self.U1 = U1
        self.W1 = W1
        self.U2 = U2
        self.W2 = W2


# section 5.2 of [Revisiting Nova on cycle of curves]()
class IVC:
    def __init__(self, F1, F2, E1, E2, func1, z1, func2, z2):
        self.F1 = F1
        self.F2 = F2
        self.E1 = E1
        self.E2 = E2
        self.R1 = [self.E1.random_element() for _ in range(10)]
        self.R2 = [self.E2.random_element() for _ in range(10)]
        self.G1 = pallas_gen
        self.G2 = vesta_gen
        self.zero_commit_1 = commit1(self.F1, self.R1, self.G1, self.F1(0))
        self.zero_commit_2 = commit2(self.F2, self.R2, self.G2, self.F2(0))

        return self.base_case(func1, z1, func2, z2)

    def base_case(self, func1, z1, func2, z2):
        """
        base case: i=0
        """
        i1 = self.F1(0)
        i2 = self.F2(0)

        # initialization
        R1_accum = CommittedRelaxedR1CS(
            self.E1(0), self.F1(0), self.E1(0), (self.F1(0), self.F1(0))
        )
        R2_accum = CommittedRelaxedR1CS(
            self.E2(0), self.F2(0), self.E2(0), (self.F2(0), self.F2(0))
        )
        R2_x0 = hash1(self.F1, (self.F1(0), z1, z1, R2_accum))
        R2_x1 = hash2(self.F2, (self.F2(0), z2, z2, R1_accum))
        R2_u0 = CommittedRelaxedR1CS(
            self.zero_commit_2, self.F2(1), self.zero_commit_2, (R2_x0, R2_x1)
        )

        # function computation
        R1_zi = func1(z1)

        # witness relation
        R1_witness_relation = (i1, z1, z1, R1_accum, R2_u0, self.zero_commit_2)
        R1_relaxed_witness = 0  # TODO: check this
        R1_w_comm = commit1(self.F1, self.R1, self.G1, R1_witness_relation)
        R2_W_iplus1 = R1_w_comm  # TODO: check this

        R1_u_Ecomm = self.zero_commit_1
        R1_ui_mu = self.F1(1)
        R1_ui_x0 = self.F1(R2_u0.x[1])  # converting to F1 since R2.x0 is in F2
        R1_ui_x1 = hash1(self.F1, (self.F1(1), z1, R1_zi, R2_accum))
        R1_ui_Wcomm = R1_w_comm
        R1_u_iplus1 = CommittedRelaxedR1CS(
            R1_u_Ecomm, R1_ui_mu, R1_ui_Wcomm, (R1_ui_x0, R1_ui_x1)
        )
        R1_w_iplus1 = (self.F1(0), R1_witness_relation)

        # function compuatation
        R2_zi = func2(z2)

        # witness relation
        R2_witness_relation = (i2, z2, z2, R2_accum, R1_u_iplus1, self.zero_commit_2)
        R2_w_comm = commit2(self.F2, self.R2, self.G2, R2_witness_relation)
        R1_accum_iplus1 = R1_u_iplus1
        R1_W_iplus1 = R1_w_iplus1

        R2_u_Ecomm = self.zero_commit_2
        R2_ui_mu = self.F2(1)
        R2_ui_x0 = self.F2(R1_u_iplus1.x[1])  # converting to F2 since R1.x0 is in F1
        R2_ui_x1 = hash2(self.F2, (self.F2(1), z2, R2_zi, R1_accum_iplus1))
        R2_u_iplus1 = CommittedRelaxedR1CS(
            R2_u_Ecomm, R2_ui_mu, R2_w_comm, (R2_ui_x0, R2_ui_x1)
        )
        R2_w_iplus1 = (self.F2(0), R2_witness_relation)

        return (
            R1_zi,
            R2_zi,
            IVCProof(
                R2_u_iplus1, R2_w_iplus1, R1_accum, R1_W_iplus1, R2_accum, R2_W_iplus1
            ),
        )

    def Prove(self, i, func1, func2, z1_0, z2_0, z1_i, z2_i, ivc_proof: IVCProof):
        # R1CS(1)
        R2_accum_iplus1 = CommittedRelaxedR1CS(
            self.E2.random_point(),
            self.F2.random_element(),
            self.E2.random_point(),
            (self.F2.random_element(), self.F2.random_element()),
        )
        R2_W_iplus1 = 0  # TODO: compute R1CS for new witness

        R2_T = 0  # TODO
        R1_witness_relation = (self.F1(i), z1_0, z1_i, ivc_proof.U2, ivc_proof.u2, R2_T)
        R1_w_comm = commit1(self.F1, self.R1, self.G1, R1_witness_relation)

        R1_z_iplus1 = func1(z1_i)

        R1_ui_x0 = self.F1(ivc_proof.u2.x[1])
        R1_ui_x1 = hash1(self.F1, (self.F1(i + 1), z1_0, R1_z_iplus1, R2_accum_iplus1))
        R1_u_iplus1 = CommittedRelaxedR1CS(
            self.zero_commit_1, self.F1(1), R1_w_comm, (R1_ui_x0, R1_ui_x1)
        )
        R1_w_iplus1 = (0, R2_W_iplus1)  # TODO: this is not correct, check this

        # R1CS(2)
        R1_accum_iplus1 = CommittedRelaxedR1CS(
            self.E1.random_point(),
            self.F1.random_element(),
            self.E1.random_point(),
            (self.F1.random_element(), self.F1.ranomd_element()),
        )
        R1_W_iplus1 = 0  # TODO: compute R1CS

        R1_T = 0  # TODO:
        R2_witness_relation = (
            self.F2(i),
            z2_0,
            z2_i,
            ivc_proof.U1,
            R1_accum_iplus1,
            R1_T,
        )
        R2_w_comm = commit2(self.F2, self.R2, self.G2, R2_witness_relation)

        R2_z_iplus1 = func2(z2_i)

        R2_ui_x0 = self.F1(R1_u_iplus1.x[1])
        R2_ui_x1 = hash2(self.F2, (self.F2(i + 1), z2_0, R2_z_iplus1, R1_accum_iplus1))
        R2_u_iplus1 = CommittedRelaxedR1CS(
            self.zero_commit_2, self.F2(1), R2_w_comm, (R2_ui_x0, R2_ui_x1)
        )
        R2_w_iplus1 = (0, R1_W_iplus1)

        return (
            R1_z_iplus1,
            R2_z_iplus1,
            IVCProof(
                R2_u_iplus1,
                R2_w_iplus1,
                R1_accum_iplus1,
                R1_W_iplus1,
                R2_accum_iplus1,
                R2_W_iplus1,
            ),
        )

    def Verify(self):
        pass


W1, X1 = witness_instance_pair(F1, 2, 3, 1)
Z1 = vector(W1 + X1)
W2, X2 = witness_instance_pair(F1, 4, 6, 0)
Z2 = vector(W2 + X2)


class TestR1CS(unittest.TestCase):
    def test_r1cs(self):
        W1, X1 = witness_instance_pair(F1, 2, 3, 1)
        Z1 = vector(W1 + X1)
        print(f"Z1: {Z1}")
        assert hadamard_prod(F1, A * Z1, B * Z1) == C * Z1

        W2, X2 = witness_instance_pair(F1, 4, 6, 0)
        Z2 = vector(W2 + X2)
        print(f"Z2: {Z2}")
        assert hadamard_prod(F1, A * Z2, B * Z2) == C * Z2


class TestNovaNIFS(unittest.TestCase):
    def test_nifs(self):
        W1, X1 = witness_instance_pair(F1, 2, 3, 1)
        Z1 = vector(W1 + X1)
        W2, X2 = witness_instance_pair(F1, 4, 6, 0)
        Z2 = vector(W2 + X2)

        nova_nifs = NIFS(F1, pallas, A, B, C)
        r = F1.random_element()
        E, W, comm_T, r1, r2 = nova_nifs.prover(Z1, Z2, r)
        (
            comm_E,
            μ,
            comm_W,
            X,
        ) = nova_nifs.verifier(comm_T, r, r1, r2)

        Z = vector(list(W) + list(X))

        assert hadamard_prod(F1, A * Z, B * Z) == μ * C * Z + E
