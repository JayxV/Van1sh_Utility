## Cry
### equivalent
```python
from Crypto.Util.number import *
from collections import namedtuple


SecretKey = namedtuple('SecretKey', ['s', 'e', 'p'])


def orthogonal_lattice(bs):

    k = bs.nrows()
    n = len(bs[0])

    _ = 1
    for each in bs:
        _ *= each.norm()

    _ = abs(_)

    g = int(2^(n - 1/2 + (n - k)*(n - k - 1)/4) * _) + 1

    M = block_matrix(ZZ, [
        Matrix(ZZ, [g * _ for _ in bs]),
        identity_matrix(n)
    ], nrows=2, subdivide=False).T

    M = M.LLL()

    return (M[:n - k].T)[-n:].T


def dec(c, sk):
    e_inv = inverse(sk.e, sk.p)
    m = (e_inv * c % sk.p) % 2
    return m


def decrypt(cip, sk):
    bits = [dec(c, sk) for c in cip]
    msg = long_to_bytes(int(''.join(map(str, bits)), 2))
    return msg


pks = []

cip = []

a = Matrix(ZZ, [pks])
a_ort = orthogonal_lattice(a)

L1 = a_ort[:-1]
L1_ort = orthogonal_lattice(L1)


u1, u2 = L1_ort[0], L1_ort[1]

while True:
    x1 = randrange(-80, 80)
    x2 = randrange(-80, 80)
    y1 = randrange(-80, 80)
    y2 = randrange(-80, 80)

    s = x1 * u1 + x2 * u2
    k = y1 * u1 + y2 * u2

    A = Matrix(ZZ, [
        s,
        k
    ])
    
    try:
        e, p = A.solve_left(a)[0]
    except:
        continue

    if gcd(p, e) == 1:
        if all([_ > 0 for _ in list(s)]) and all([(_ % 2) == 1 for _ in list(s)]) and p > sum(list(s)):
            sk = SecretKey(list(s), int(e), int(p))
            print(sk)
            print(decrypt(cip, sk))
            raise Exception('GG')

```
### d3bug
给了70个比特 但是未知量只有64个，列70个64元GF(2)上的方程解方程即可。
```python
flag = var(["f%d" % i for i in range(64)])
variables = list(flags)
mask = Integer(0b1010010000001000000010001001010010100100000010000000100010010100).digits(2, padto=64)[::-1]
standard = "01111101111010111000010010111001101"
my = "00100110001000110001101010101001001"

for i in range(35):
    variables.append(sum([each1 * each2 for each1, each2 in zip(variables[-64:], mask)]))

L1 = variables[-35:]
variables = list(flags)
for i in range(35):
    variables.append(sum(variables[-64:]))

L2 = variables[-35:]
y = []
A = []
for i in range(35):
    this1 = []
    this2 = []
    for j in range(64):
        this1.append(int(L1[i].coefficient(flags[j])) % 2)
        this2.append(int(L2[i].coefficient(flags[j])) % 2)
    A.append(this1)
    A.append(this2)
    y.append(int(standard[i]))
    y.append(int(my[i]))

A = Matrix(GF(2), A)
y = vector(GF(2), y)

res = A.solve_right(y)

from Crypto.Util.number import long_to_bytes
print(long_to_bytes(int(''.join([str(each) for each in res]), 2)))
```

### d3qcg
二元coppersmith
```python
import itertools

def small_roots(f, bounds, m=1, d=None):
	if not d:
		d = f.degree()

	R = f.base_ring()
	N = R.cardinality()
	
	f /= f.coefficients().pop(0)
	f = f.change_ring(ZZ)

	G = Sequence([], f.parent())
	for i in range(m+1):
		base = N^(m-i) * f^i
		for shifts in itertools.product(range(d), repeat=f.nvariables()):
			g = base * prod(map(power, f.variables(), shifts))
			G.append(g)

	B, monomials = G.coefficient_matrix()
	monomials = vector(monomials)

	factors = [monomial(*bounds) for monomial in monomials]
	for i, factor in enumerate(factors):
		B.rescale_col(i, factor)

	B = B.dense_matrix().LLL()

	B = B.change_ring(QQ)
	for i, factor in enumerate(factors):
		B.rescale_col(i, 1/factor)

	H = Sequence([], f.parent().change_ring(QQ))
	for h in filter(None, B*monomials):
		H.append(h)
		I = H.ideal()
		if I.dimension() == -1:
			H.pop()
		elif I.dimension() == 0:
			roots = []
			for root in I.variety(ring=ZZ):
				root = tuple(R(root[var]) for var in f.variables())
				roots.append(root)
			return roots

	return []

key = {'a': 3591518680290719943596137190796366296374484536382380061852237064647969442581391967815457547858969187198898670115651116598727939742165753798804458359397101, 'c': 6996824752943994631802515921125382520044917095172009220000813718617441355767447428067985103926211738826304567400243131010272198095205381950589038817395833, 'p': 7386537185240346459857715381835501419533088465984777861268951891482072249822526223542514664598394978163933836402581547418821954407062640385756448408431347}
a, c, p = key.values()

hint = [67523583999102391286646648674827012089888650576715333147417362919706349137337570430286202361838682309142789833, 70007105679729967877791601360700732661124470473944792680253826569739619391572400148455527621676313801799318422]
Bits = 512
UnKnownBits = 146

PR.<x, y> = PolynomialRing(Zmod(p))
f = (a * ((hint[0] << UnKnownBits) + y)^2 + c) - ((hint[1] << UnKnownBits) + x)
bounds = (2^(UnKnownBits + 1), ) * 2

x, y = small_roots(f, bounds, m=4, d=4)[0]
y = (hint[0] << UnKnownBits) + y
data = (y - c) * inverse_mod(a, p) % p
num = pow(data, (p + 1) // 4, p)
enc = 6176615302812247165125832378994890837952704874849571780971393318502417187945089718911116370840334873574762045429920150244413817389304969294624001945527125
from hashlib import sha512
from Crypto.Util.number import *
flag = bytes_to_long(sha512(b'%d'%(num)).digest())^^enc
print(long_to_bytes(flag))
```
### d3factor

```python
from Crypto.Util.number import *
from gmpy2 import iroot
from hashlib import md5


N = 
e1 = 
e2 = 

e1, e2 = e2 ,e1

PR.<x> = Zmod(N)[]

f = e1 * e2 * x - (e2 - e1)
f = f.monic()

res = f.small_roots(X=2^1000, beta=0.125, epsilon=0.02)[0]

p = GCD(Integer(f(res)), N)
p = iroot(p, 6)[0]

q = N // p^7

assert p > q

e = 0x10001
phi = (p - 1) * (q - 1)

c = 
d = inverse(e, phi)
msg = long_to_bytes(pow(c, d, p * q))

print(msg)

flag = 'd3ctf{' + md5(msg).hexdigest() + '}'

print(flag)

```
