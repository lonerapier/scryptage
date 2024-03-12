# adapted from:
# - https://montekki.github.io/thaler-ch4-4/
# - https://github.com/iammadab/zk/blob/main/protocols/src/gkr

from sage.all import *
from sage.rings.polynomial.polydict import PolyDict
import unittest

load("sumcheck.sage")
load("poly_utils.sage")

p = 21888242871839275222246405745257275088696311157297823662689037894645226208583
F = GF(p)

add_gate = "add"
mul_gate = "mul"


def lagrange_basis(bin_i, num):
    res = multiplicative_identity(num)
    zero_tuple = [0] * num

    for j in range(num):
        indeterminate = zero_tuple.copy()
        indeterminate[j] = 1
        if bin_i[j] == "0":
            # print("res:", res)
            # print("poly:", PolyDict({tuple(zero_tuple): 1, tuple(indeterminate): -1}))

            res *= PolyDict({tuple(zero_tuple): 1, tuple(indeterminate): -1})
        else:
            res *= PolyDict({tuple(indeterminate): 1})

    return res


def mle(fw):
    num = log(len(fw), 2)
    if num == 0:
        num = 1
    delta_w = additive_identity(num)
    # print("delta_w:", delta_w)
    # print("num:", num, "len(fw):", len(fw))
    for i in range(len(fw)):
        # print("delta_w:", delta_w)
        # print(
        #     "lagrange_basis:", lagrange_basis(to_bin(i, num), num).scalar_lmult(fw[i])
        # )
        delta_w += lagrange_basis(to_bin(i, num), num).scalar_lmult(fw[i])

    return delta_w


class Gate:
    """
    Gate class to represent a fan-in 2 gate for the circuit
    """

    def __init__(self, gate_type, output, input_1, input_2):
        self.type = gate_type
        self.input_1 = input_1
        self.input_2 = input_2
        self.output = output

    def to_bin(self, layer_index):
        return (
            to_bin(self.output, layer_index)
            + to_bin(self.input_1, layer_index + 1)
            + to_bin(self.input_2, layer_index + 1)
        )


class Layer:
    def __init__(self, add_gates, mul_gates, index):
        self.add_gates = add_gates
        self.mul_gates = mul_gates
        self.index = index

    def len(self):
        return len(self.add_gates) + len(self.mul_gates)

    def max_gate_index(self):
        return max(
            [max(gate.input_1, gate.input_2, gate.output) for gate in self.add_gates]
            + [max(gate.input_1, gate.input_2, gate.output) for gate in self.mul_gates]
        )

    def add_mul_mle(self, layer_len):
        """
        multilinear extension polynomial for add and mul gates in a layer.
        # inputs:
        - `layer_len`: total gates in layer

        TODO: add support for non-uniform circuits
        """

        if len(self.add_gates) == 0 and len(self.mul_gates) == 0:
            return additive_identity(1), additive_identity(1)

        layer_num = log(layer_len, 2)

        add_mle = additive_identity(layer_num + (2 * layer_num + 1))
        mul_mle = additive_identity(layer_num + (2 * layer_num + 1))

        for gate in self.add_gates:
            gate_bin_string = gate.to_bin(
                layer_num,
            )
            add_mle += lagrange_basis(gate_bin_string, len(gate_bin_string))

        for gate in self.mul_gates:
            gate_bin_string = gate.to_bin(layer_num)
            mul_mle += lagrange_basis(gate_bin_string, len(gate_bin_string))

        return add_mle, mul_mle


class Circuit:
    def __init__(self, layers):
        self.layers = layers

    def depth(self):
        return len(self.layers)

    def add_mul_mle(self, layer_index):
        """
        multilinear extension polynomial for add and mul gates in a layer.
        assumption: circuit is uniform: i.e. next layer len is 2*current_layer_len

        inputs:
        - `layer_index`: index of the layer

        outputs:
        - `add_mle`: multilinear extension polynomial for add gates
        - `mul_mle`: multilinear extension polynomial for mul gates
        """

        return self.layers[layer_index].add_mul_mle(self.layers[layer_index].len())

    def evaluate(self, inputs):
        """
        evaluate the circuit in reverse and return the evaluations layer wise
        """

        assert len(inputs) == self.layers[len(self.layers) - 1].max_gate_index() + 1

        evaluations = [inputs]
        # go through each layer in reverse and evaluate add and mul gate to find evaluations
        for i in range(len(self.layers) - 1, -1, -1):
            layer_evaluations = [F(0)] * pow(
                2, i
            )  # TODO: works for uniform but not for non-uniform circuits
            for gate in self.layers[i].add_gates:
                layer_evaluations[gate.output] = (
                    inputs[gate.input_1] + inputs[gate.input_2]
                )

            for gate in self.layers[i].mul_gates:
                layer_evaluations[gate.output] = (
                    inputs[gate.input_1] * inputs[gate.input_2]
                )

            inputs = layer_evaluations

            evaluations.append(layer_evaluations)

        evaluations.reverse()

        return evaluations

    def w(self, evaluations, layer_index):
        """
        multilinear extension of w polynomial

        Inputs:
        - `evaluations`: list of evaluations of the circuit at layer level
        - `layer_index`: index of the layer

        Outputs:
        - `w_mle`: multilinear extension polynomial for w
        """
        assert len(evaluations) == len(self.layers) + 1
        assert layer_index <= len(self.layers)
        return mle(evaluations[layer_index])


class GateEvalExtension:
    def __init__(self, add_mle, mul_mle, w_b, w_c, randomness):
        self.add_mle = add_mle
        self.mul_mle = mul_mle
        self.w_b = w_b
        self.w_c = w_c
        self.randomness = randomness

    def sumcheck_prove(self):
        pass


def q_function(line_func, G, w):
    """
    create a univariate polynomial equal to w(l(x)).
    substitutes `l[i]` for each ith variable in multivariate poly `w`
    """
    def evaluate_term_with_poly(term, G, lines):
        poly = G[0]
        for i in range(len(term)):
            poly = poly + lines[i] ** term[i]
        return poly

    terms = w.dict()
    assert len(line_func) == num_vars(w)
    unipoly = G[0]
    for term, coeff in terms.items():
        unipoly_term = evaluate_term_with_poly(term, G, line_func)
        unipoly = unipoly + (coeff*unipoly_term)

    return unipoly

def verify_sumcheck(proof, r, sum):
    """
    verify sumcheck non-interactively leaving final sumcheck sum proof
    """
    for i in range(len(proof)):
        unipoly = proof[i]
        expected = unipoly(0) + unipoly(1)
        assert(expected == sum)
        # assert(unipoly.degree() <= lookup_degree_table[i])

        sum = unipoly(r[i])
    return r, sum

class GKR:
    def __init__(self, F, circuit: Circuit, inputs, evaluations):
        """
        inputs:
        - `circuit`: circuit object
        - `evaluations`: list of evaluations of the circuit at layer level
        """
        self.F = F
        G.<x> = F['x'] # unipoly
        self.G = G
        self.x = x # poly invariate
        self.circuit = circuit
        self.inputs = inputs
        self.evaluations = evaluations
        self.randomness = []
        self.round_sum = None

    def line_func(self, b, c):
        """
        generate a line such that l(0) = b, l(1) = c
        Note: b and c are vectors, so line func is also a list of univariate polys
        """

        assert len(b) == len(c)
        lines = []
        for i in range(len(b)):
            unipoly = b[i] + self.x * (c[i]-b[i])
            lines.append(unipoly)
        return lines

    def evaluate_line_func(self, lines, r):
        assert len(lines) == len(r)
        result = []
        for i in range(len(r)):
            value = lines[i](r[i])
            result.append(value)
        return result

    def prove_output_mle(self):
        w_0 = self.circuit.w(self.evaluations, 0)
        # m = evaluate(w_0, self.F, randomness)

        # self.randomness = randomness

        return w_0

    def prove(self, layer_index, randomness_i1):
        """
        prove the GKR protocol interactively. Prover will prove the sumcheck protocol for each layer and


        inputs:
        - `layer_index`: index of the layer being proven
        - `randomness`: random field elements for w_i sumcheck

        outputs:
        - `sumcheck_proofs`: non-interactive sumcheck proofs (unipolys)
        - `sumcheck_proof_sum`: final round sumcheck proof sum
        - `sumcheck_randomness`: random variables used in sumcheck protocol. b* and c* of round polynomial
        - `q_function`: w(l(x)), a univariate polynomial of degree k_{i+1} which restrict w poly to line l
        """

        # if layer_index == None:
        #     return None, None, None, None, self.circuit.w(self.evaluations, 0)
        add_mle_i, mul_mle_i = self.circuit.add_mul_mle(layer_index)
        w_i1 = self.circuit.w(self.evaluations, layer_index+1)

        # round_poly_poly = (add_mle * (w+w)) + (mul_mle * (w*w))
        round_poly = GateEvalExtension(add_mle_i, mul_mle_i, w_i1, w_i1, self.randomness)

        sumcheck_proof, sumcheck_randomness = round_poly.sumcheck_prove() # TODO: sumcheck_randomness should come from verifier

        b_star, c_star = (
                sumcheck_randomness[: len(sumcheck_randomness) / 2],
                sumcheck_randomness[len(sumcheck_randomness) / 2 :],
            )
        line_function = self.line_func(b_star, c_star) # list of univariate line polys

        q = q_function(line_function, self.G, w_i1)
        m = evaluate(w_i1, self.F, randomness_i1)
        self.randomness = randomness_i1

        return sumcheck_proof, m, sumcheck_randomness, q

    def verify(self):
        """
        verify the GKR protocol
        """

        output_mle = self.prove_output_mle() # w_0
        r_i = [self.F.random_element() for i in range(n_vars(output_mle))] # r_0
        self.randomness = r_i
        m = evaluate(output_mle, self.F, r_i) # m = D(r_0)

        for i in range(self.circuit.depth()):
            # prove round `i`
            (
                sumcheck_proofs,
                sumcheck_proof_sum,
                sumcheck_randomness, # b* and c*
                q_function,
            ) = self.prove(i, r_i)

            assert m == sumcheck_proof_sum
            assert len(sumcheck_randomness) % 2 == 0

            # verify sumcheck proof
            add_mle, mul_mle = self.circuit.add_mul_mle(i)
            sumcheck_randomness, verify_sum = verify_sumcheck(sumcheck_proofs, sumcheck_randomness, sumcheck_proof_sum)

            # final check of sumcheck protocol
            #
            # evaluate f_ri_b_c = add(r, b*, c*) * (w(b*) + w(c*)) + mul(r, b*, c*) * (w(b*) * w(c*))
            # at randomness sent by prover
            # q = w(l(x))
            # w(b) = q(0)
            # w(c) = q(1)
            w_b = q_function(0)
            w_c = q_function(1)
            evaluation_input = r_i
            evaluation_input.extend(sumcheck_randomness)
            add_result = evaluate(add_mle, evaluation_input) * (w_b + w_c)
            mul_result = evaluate(mul_mle, evaluation_input) * (w_b * w_c)
            round_poly_eval = add_result + mul_result

            # assert final sumcheck sum
            assert round_poly_eval == verify_sum

            # reduce two points to one line, get line function, r_i+1, m_i+1
            b, c = (
                sumcheck_randomness[: len(sumcheck_randomness) / 2],
                sumcheck_randomness[len(sumcheck_randomness) / 2 :],
            )
            l_i1 = self.line_func(b, c)
            r_star = self.F.random_element()
            r_i = self.evaluate_line_func(l_i1, r_star)
            m = q_function(r_star) # prover will prove this value in next round


class TestMLE(unittest.TestCase):
    def test_mle(self):
        print("starting mle test")
        fw = [1, 2, 8, 10]
        assert evaluate(mle(fw), F, [1, 1]) == 10

        fw = [1, 2, 8, 10, 2, 4, 16, 20]
        assert evaluate(mle(fw), F, [1, 1, 1]) == 20

    def test_lagrange_basis(self):
        print("starting lagrange basis test")
        assert lagrange_basis(to_bin(0, 2), 2) == PolyDict(
            {(0, 0): 1, (0, 1): -1, (1, 0): -1, (1, 1): 1}
        )
        assert lagrange_basis(to_bin(3, 2), 2) == PolyDict({(1, 1): 1})


def sample_uniform_circuit():
    layer_0 = Layer([], [Gate(mul_gate, 0, 0, 1)], 0)
    layer_1 = Layer([Gate(add_gate, 0, 0, 1)], [Gate(mul_gate, 1, 2, 3)], 1)
    layer_2 = Layer(
        [Gate(add_gate, 0, 0, 1), Gate(add_gate, 1, 2, 3), Gate(add_gate, 2, 4, 5)],
        [Gate(mul_gate, 3, 6, 7)],
        2,
    )
    return Circuit([layer_0, layer_1, layer_2])


class TestCircuit(unittest.TestCase):
    def test_gate(self):
        print("starting gate test")
        gate = Gate(add_gate, 0, 0, 1)
        assert gate.to_bin(1) == ["0", "0", "0", "0", "1"]

        gate = Gate(mul_gate, 1, 2, 3)
        assert gate.to_bin(1) == ["1", "1", "0", "1", "1"]

    def test_layer(self):
        print("starting layer test")
        layer = Layer([Gate(add_gate, 0, 0, 1)], [Gate(mul_gate, 1, 2, 3)], 1)
        add_mle, mul_mle = layer.add_mul_mle(2)

        assert evaluate(mul_mle, F, [1, 1, 0, 1, 1]) == 1
        assert evaluate(add_mle, F, [0, 0, 0, 0, 1]) == 1

    def test_circuit(self):
        print("starting circuit test")
        sample_circuit = sample_uniform_circuit()
        assert sample_circuit.depth() == 3

        add_mle_0, mul_mle_0 = sample_circuit.add_mul_mle(0)
        assert evaluate(mul_mle_0, F, [0, 0, 1]) == 1
        assert evaluate(add_mle_0, F, [0]) == 0

        add_mle_1, mul_mle_1 = sample_circuit.add_mul_mle(1)
        assert evaluate(add_mle_1, F, [0, 0, 0, 0, 1]) == 1
        assert evaluate(mul_mle_1, F, [1, 1, 0, 1, 1]) == 1

        add_mle_2, mul_mle_2 = sample_circuit.add_mul_mle(2)
        assert evaluate(add_mle_2, F, [0, 0, 0, 0, 0, 0, 0, 1]) == 1
        assert evaluate(mul_mle_2, F, [1, 1, 1, 1, 0, 1, 1, 1]) == 1

        evaluations = sample_circuit.evaluate([1, 2, 3, 4, 5, 6, 7, 8])
        assert evaluations[0][0] == 6160
        assert evaluations[1][0] == 10
        assert evaluations[1][1] == 616
        assert evaluations[2][0] == 3
        assert evaluations[2][1] == 7

        w_0 = sample_circuit.w(evaluations, 0)
        assert evaluate(w_0, F, [0]) == 6160

        w_1 = sample_circuit.w(evaluations, 1)
        assert evaluate(w_1, F, [0]) == 10
        assert evaluate(w_1, F, [1]) == 616

        w_2 = sample_circuit.w(evaluations, 2)
        assert evaluate(w_2, F, [0, 0]) == 3

def sample_uniform_gkr():
    sample_circuit = sample_uniform_circuit()
    input = [F(1), F(2), F(3), F(4), F(5), F(6), F(7), F(8)]
    evaluations = sample_circuit.evaluate(input)

    gkr = GKR(F, sample_circuit, input, evaluations)

    return gkr, sample_circuit, input, evaluations

class TestGKR(unittest.TestCase):
    def test_line_func(self):
        gkr, circuit, inputs, evaluations = sample_uniform_gkr()

        num = 10
        b = [F.random_element() for i in range(num)]
        c = [F.random_element() for i in range(num)]

        line_func = gkr.line_func(b, c)
        assert line_func[0](0) == b[0]
        assert line_func[0](1) == c[0]


    def test_gkr(self):
        print("starting gkr test")

        sample_circuit = sample_uniform_circuit()
        input = [F(1), F(2), F(3), F(4), F(5), F(6), F(7), F(8)]
        evaluations = sample_circuit.evaluate(input)

        gkr = GKR(F, sample_circuit, input, evaluations)

        # gkr.prove()

        # gkr.verify()


if __name__ == "__main__":
    unittest.main()
