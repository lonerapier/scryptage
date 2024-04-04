# adapted from [Anatomy of a STARK](https://aszepieniec.github.io/stark-anatomy/fri)
# and [winterfell FRI](https://github.com/facebook/winterfell/blob/main/fri)

from sage.all import *
import unittest
from hashlib import blake2b


def evaluate_univ_poly(poly, domain):
    return [poly(domain[i]) for i in range(len(domain))]


class Merkle:
    def __init__(self, evaluations):
        self.leaves = evaluations
        self.num_layers = log2(evaluations)
        self.tree = self.build_tree(evaluations)

    def build_tree(self, evlautions):
        return [[]]

    def root(self):
        return self.tree[0][0]

    def open(self, leaf_index):
        pass

    def verify(root, path, leaf):
        pass


class FRILayer:
    def __init__(self, tree: Merkle, evaluations) -> None:
        self.tree = tree
        self.evaluations = evaluations


class FRI:
    def __init__(
        self, F, blowup_factor, omega, offset, domain_length, num_queries
    ) -> None:
        self.F = F
        self.blowup_factor = blowup_factor
        self.omega = omega
        self.offset = offset
        self.domain_length = domain_length
        self.num_rounds = self.calc_rounds(domain_length)
        self.num_queries = num_queries

    def calc_rounds(self, domain_length):
        codeword_length = domain_length
        num_rounds = 0
        while (
            codeword_length > self.blowup_factor
            and 4 * self.num_queries < codeword_length
        ):
            num_rounds += 1
            codeword_length >>= 1
        return num_rounds

    def sample_index(self, digest, size):
        pass

    def sample_indices(self, size, reduced_size, number):
        assert (
            number <= 2 * reduced_size
        ), "not enough entropy in indices wrt last codeword"
        assert (
            number <= reduced_size
        ), f"cannot sample more indices than available in last codeword; requested: {number}, available: {reduced_size}"

        indices = []
        reduced_indices = []
        counter = 0
        while len(indices) < number:
            index = self.sample_index(blake2b(seed + bytes(counter)).digest(), size)
            reduced_index = index % reduced_size
            counter += 1
            if reduced_index not in reduced_indices:
                indices += [index]
                reduced_indices += [reduced_index]

        return indices

    def commit(self, codeword, transcript):
        one = self.F(1)
        two_inv = 1 / self.F(2)
        omega = self.omega
        offset = self.offset

        codewords = []
        layers = []
        alphas = []

        for i in range(self.num_rounds):
            codewords.append(codeword)

            # merkelise codeword and add root
            evaluation_tree = Merkle(codeword)
            layers.append(FRILayer(tree=evaluation_tree, evaluations=codeword))

            root = evaluation_tree.root()
            transcript.append(root)

            # don't apply drp after last round
            if i == (self.num_rounds - 1):
                break

            # apply DRP
            alpha = self.F.random_element()
            alphas.append(alpha)

            codeword = [
                two_inv
                * (
                    (one + alpha / (offset * (omega ^ k))) * codeword[k]
                    + (one - alpha / (offset * (omega ^ k)))
                    * codeword[len(codeword) // 2 + k]
                )
                for k in range(len(codeword) // 2)
            ]

            omega <<= 1
            offset <<= 1

        # take remainder codeword for polynomial interpolation check
        transcript.append(codeword)
        codewords.append(codeword)

        return codewords, layers, alphas

    def query(self, transcript, indices, layer_prev, layer_curr):
        codeword_prev = layer_prev.evaluations
        codeword = layer_curr.evaluations

        a_indices = [index for index in indices]
        b_indices = [index + len(codeword_prev) // 2 for index in a_indices]

        for i in range(self.num_queries):
            transcript.append(
                (
                    codeword_prev[a_indices[i]],
                    codeword_prev[b_indices[i]],
                    codeword[indices[i]],
                )
            )
            transcript.append(layer_prev.tree.open(a_indices[i]))
            transcript.append(layer_prev.tree.open(b_indices[i]))
            transcript.append(layer_curr.tree.open(indices[i]))

    def Prove(self, codeword):
        transcript = []
        codewords, layers, alphas = self.commit(codeword, transcript)

        colinearity_indices = self.sample_indices(
            len(codewords[1]), len(codewords[-1]), self.num_queries
        )

        for i in range(len(codewords) - 1):
            indices = [index % len(codewords[i]) // 2 for index in colinearity_indices]
            self.query(transcript, indices, layers[i], layers[i + 1])

        return transcript, alphas, colinearity_indices

    def test_colinearity(self, p1, p2, p3):
        R = PolynomialRing(self.F, "x")
        poly = R.lagrange_polynomial([p1, p2])
        return poly(p3[0]) == p3[1]

    def Verify(self, transcript, alphas, colinearity_indices):
        """
        verifes an FRI proof that codeword is close to a RS codeword.
        Outputs `true/false`.

        input:
        - `transcript`: list containing proof elements like roots of merkle tree, paths to query indices
        - `alphas`: randomness supplied by verifer to fold codewords
        - `colinearity_indices`: indices for query phase of verifier
        """
        omega = self.omega
        offset = self.offset

        # get all roots sent by prover
        roots = [transcript.pop(0) for _ in range(self.num_rounds)]

        last_codeword = transcript.pop(0)

        # check if last codeword matches with root sent by prover
        root = Merkle(last_codeword)
        assert root == roots[-1]

        last_omega = omega << (self.num_rounds - 1)
        last_offset = offset << (self.num_rounds - 1)

        # check if last codeword is low degree
        last_codeword_degree = (len(last_codeword) // self.blowup_factor) - 1
        R = PolynomialRing(self.F, "x")

        last_poly_domain = [
            (last_offset * (last_omega ^ i)) for i in range(len(last_codeword))
        ]
        last_codeword_domain_and_points = [
            (last_offset * (last_omega ^ i), last_codeword[i])
            for i in range(len(last_codeword))
        ]
        poly = R.lagrange_polynomial(last_codeword_domain_and_points)
        poly_evaluation = evaluate_univ_poly(poly, last_poly_domain)
        assert poly_evaluation == last_codeword
        assert poly.degree() < last_codeword_degree

        # start query round, for each folding step, do queries and check at random points
        initial_codeword_length = len(last_codeword) << ((self.num_rounds) - 1)
        for k in range(self.num_rounds - 1):
            indices = [
                index[i] % (initial_codeword_length >> (k + 1))
                for index in colinearity_indices
            ]
            a_indices = [index for index in indices]
            b_indices = [
                index + (initial_codeword_length >> (k + 1)) for index in a_indices
            ]

            # do colinearity test and check merkle path for each query
            for j in self.num_queries:
                (ay, by, cy) = transcript.pop(0)
                ax = offset * (omega ^ a_indices[j])
                bx = offset * (omega ^ b_indices[j])
                cx = alphas[k]

                assert self.test_colinearity((ax, ay), (bx, by), (cx, cy))

                a_path = transcript.pop(0)
                assert Merkle.verify(roots[k], a_path, ay)

                b_path = transcript.pop(0)
                assert Merkle.verify(roots[k], b_path, by)

                c_path = transcript.pop(0)
                assert Merkle.verify(roots[k + 1], c_path, cy)

            omega <<= 1
            offset <<= 1
