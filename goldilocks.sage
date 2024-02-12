# https://twitter.com/cronokirby/status/1565381285169741829
# https://xn--2-umb.com/22/goldilocks/

p = 2 ^ 64 - 2 ^ 32 + 1
print(p)
F = GF(p)

print(F.primitive_element())

g = F(7)
print(g.multiplicative_order().factor())

rou = 2 ^ 32
rou2 = 2 ^ 16 + 1
rou3 = 2 ^ 8 + 1
p2 = g ^ int((p - 1) / rou)
g2 = F(p2)  # gen for 2^32 rou
print(g2, g2.multiplicative_order().factor())


g3 = F(20033703337)
print(g3, g3.multiplicative_order().factor())

# for i in range(rou):
#     if gcd(i, rou) != 1:
#         continue
#     num = rou / i
#     tf = F(num)
#     print(i, tf.order())
