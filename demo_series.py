"""
demo_series.py

Muestra como se ven las series de Taylor (exp, sin, cos) aplicadas a
operadores de OpenKet, en distintos casos: un operador generico que no
colapsa nada, un proyector idempotente, un operador nilpotente, un
operador diagonal con autovalores distintos, y una rotacion de spin.

Correr con: python demo_series.py
"""

from sympy import var, E, cos

from openket_V2 import Ket, Bra, Operator
from Series_v2 import qexp, qsin, qcos


def titulo(texto):
    print("\n" + "=" * 70)
    print(texto)
    print("=" * 70)


# ---------------------------------------------------------------------
titulo("1. Operador generico -- la serie NO colapsa, se ve termino a termino")
# ---------------------------------------------------------------------
A = Operator("A")
for orden in range(0, 6):
    print(f"  orden {orden}:  e^A =", qexp(A, order=orden))


# ---------------------------------------------------------------------
titulo("2. Proyector idempotente |0><0|  (P^2 = P)")
# ---------------------------------------------------------------------
P = Ket(0) * Bra(0)
print("  P**2 =", P ** 2, "  (sigue siendo P, por eso la serie colapsa a un solo termino)")
print()
for orden in [2, 3, 8, 15]:
    print(f"  orden {orden:2d}:  e^P =", qexp(P, order=orden))
print()
print("  Forma cerrada conocida: e^P = 1 + (e-1)*P")
print("  e - 1 =", float(E - 1))
print(Ket(0)*qexp(P, order=1)*Bra(0))


# ---------------------------------------------------------------------
titulo("3. Operador nilpotente N = |1><0|  (N^2 = 0)")
# ---------------------------------------------------------------------
N_op = Ket(1) * Bra(0)
print("  N**2 =", N_op ** 2)
print("  e^N (orden 1)  =", qexp(N_op, order=1))
print("  e^N (orden 10) =", qexp(N_op, order=10), "  (exacto en cualquier orden >= 1)")


# ---------------------------------------------------------------------
titulo("4. Operador diagonal con autovalores DISTINTOS: N = |0><0| + 2|1><1|")
# ---------------------------------------------------------------------
Nop = Ket(0) * Bra(0) + 2 * Ket(1) * Bra(1)
print("  N =", Nop)
print("  N**2 =", Nop ** 2, "  (los cruzados |0><1|, |1><0| se anulan por ortogonalidad)")
print()
print("  La serie creciendo orden por orden:")
for orden in range(1, 7):
    print(f"    orden {orden}:", qexp(Nop, order=orden))
print()
print("  Aplicada directamente a cada ket base (asi se evita tener que sumar")
print("  el termino escalar '1' a mano):")
r30 = qexp(Nop, order=30)
resultado_0 = r30 * Ket(0)
resultado_1 = r30 * Ket(1)
print("    e^N * |0> =", resultado_0, "  ~", float(resultado_0.coef))
print("    e^N * |1> =", resultado_1, "  ~", float(resultado_1.coef))
print("  valores exactos esperados: e^1 =", float(E), "   e^2 =", float(E ** 2))


# ---------------------------------------------------------------------
titulo("5. Rotacion de spin: e^(i*theta*X) sobre la compuerta Pauli-X")
# ---------------------------------------------------------------------
from gates import X
from sympy import I

theta = var("theta")
Xgate = X()
print("  X =", Xgate)
print("  X**2 =", Xgate ** 2, "  (identidad en el subespacio de 2 niveles)")
print()
estado = qexp(I * theta * Xgate, order=8) * Ket(0)
print("  e^(i*theta*X) * |0>  (orden 8) =")
print(" ", estado)
print()
print("  Comparar el coeficiente de |0> contra la serie de cos(theta):")
print("   ", cos(theta).series(theta, 0, 9).removeO())


# ---------------------------------------------------------------------
titulo("6. qsin / qcos (mismo mecanismo, otra funcion de sympy)")
# ---------------------------------------------------------------------
print("  sin(P) (orden 10) =", qsin(P, order=10))
print("  cos(P) (orden 10) =", qcos(P, order=10))

