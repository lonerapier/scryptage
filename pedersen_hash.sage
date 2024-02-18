# 4-bit lookup window pedersen hash
# specifications from https://github.com/iden3/iden3-docs/blob/master/source/iden3_repos/research/publications/zkproof-standards-workshop-2/pedersen-hash/Pedersen-Hash.pdf
# modified from https://github.com/henryfour/pedersen-go

from math import ceil
import random


load("baby_jubjub.sage")

genPoint = FindGenPoint(Fr, EC, A, N)
lookup_bits = 4


def encode(window, lookup_bits):
    num = 1 + window[0] + 2 * window[1] + 4 * window[2]
    if window[lookup_bits - 1]:
        num = -1 * num
    return num


def pointFromBytes(EC, name, i):
    current = EC.random_point()
    try:
        if current.order != N:
            raise NotInSubgroup("point not in prime-ordered subgroup")
    except NotInSubgroup:
        current = random.randint(1, 1e3) * FindGenPoint(Fr, EC, A, N)
    return current


# def pointFromBytes():
def PedersenHash(msg, lookup_bits):
    total_windows = ceil(len(msg) / lookup_bits)
    windows = []
    i = 0
    for i in range(total_windows - 1):
        windows.append(msg[i * lookup_bits : (i + 1) * lookup_bits])
    windows.append(msg[i * lookup_bits :])

    result = NewPoint(EC)
    current = NewPoint(EC)
    for j in range(total_windows):
        if j % 50 == 0:
            current = pointFromBytes(EC, "", j / 50)

        num = encode(windows[j], lookup_bits=lookup_bits)
        segment = num * current
        result = result + segment

        current = (2 ^ (lookup_bits + 1)) * current  # recheck
