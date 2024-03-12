n = 2**55 - 16
print(n.str(16))
a = 36028797018963664
b = 36028797018963952

print(a.str(16), b.str(16))

P = (2**255) - 19
print(P.str(16))

P4 = (P << 4) % (2**256)
print(P4.str(16))

F = GF(P)
ec = EllipticCurve(F, [0, 486662, 0, 1, 0])
