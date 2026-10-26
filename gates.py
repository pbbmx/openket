from sympy import I
from openket_V2 import Ket, Bra

__all__ = ['X', 'Y', 'Z']

def X(A='default'):
    return Ket(0, A)*Bra(1, A) + Ket(1, A)*Bra(0, A)

def Y(A='default'):
    return -I*Ket(0, A)*Bra(1, A) + I*Ket(1, A)*Bra(0, A)

def Z(A='default'):
    return Ket(0, A)*Bra(0, A) - Ket(1, A)*Bra(1, A)