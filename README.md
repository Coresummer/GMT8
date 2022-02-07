Pairing Mapping Implementation with 

GMT8-542 parameter

u = 0x7452

2NAF(u) = 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1

ht = -1
hy = 0xffbbffffffffffffc020

2NAF(hy) = 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, -1, 0, 0, 1

p(u, ht ,hy) =

0x347111bfc75e57d130de7be68437
c8d75455d209459d421455023bee14
df9fe75aa4734686ca3d08c1fa5941
00d79421d56c53899ee0f066fad9eb
45b0985dbdbba2dcc1

r(u,ht,hy) = 
0xff005ff010fbfd093a41afce5a02
6f5b55729902a9b26b3179c1806080
400011

t(u,ht,hy) =
-0xff005ff010fbfd093a41afce5a0
26f5b55729902a9b26b307a0180607
c3e000e


Towering Extension Field

Fp2 : Fp[ω]/ (ω + ω^p + 1)

Fp4 : Fp2 [v]/(v^2 − (1, −1))

Fp8 : Fp4 [g]/(g^2 − v),
where (1, −1) ∈ Fp2
ω + ω^p = 1

With Elliptic Curve and Twisted Curve:

E : y^2 = x^3 + x

E': y^2 = x^3 + xv^2