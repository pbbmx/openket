"""
Series de Taylor (de sympy) aplicadas a expresiones de OpenKet.

Idea: sympy ya sabe calcular los coeficientes de Taylor de exp(x), sin(x),
etc. alrededor de x=0. Tomamos esos coeficientes (numeros/escalares) y los
usamos para construir sum_k coef_k * A**k, donde A es un operador/DiracMul/
DiracAdd de OpenKet y A**k ya usa nuestro propio __pow__ (ya validado:
recorre _evaluar como cualquier otra multiplicacion).

No reinventamos el desarrollo de Taylor -- se lo pedimos a sympy.series()
y solo hacemos la sustitucion x**k -> A**k al final.
"""

from sympy import Poly, Symbol, S, exp as sympy_exp, sin as sympy_sin, cos as sympy_cos

from openket_V2 import DiracObject

__all__ = ["qseries", "qexp", "qsin", "qcos"]


def qseries(func_expr, A, x, order):
    """
    Desarrolla func_expr(x) en serie de Taylor alrededor de x=0 hasta orden
    `order` (usando sympy.series), y sustituye cada x**k por A**k.

    :param func_expr: expresión simbólica de sympy en la variable `x`
                       (ej. exp(x), sin(x), 1/(1-x)).
    :param A: operador/expresión de OpenKet (DiracMul, DiracAdd, Operator...)
              que reemplaza a x.
    :param x: el símbolo de sympy usado en func_expr.
    :param order: número de términos a incluir (orden de truncamiento).
    :return: la expresión resultante en el álgebra de OpenKet (o un escalar
             puro de sympy si A termina siendo puramente numérico).
    """
    taylor = func_expr.series(x, 0, order + 1).removeO()
    poly = Poly(taylor, x)
    coeffs = poly.all_coeffs()[::-1]  # índice k = coeficiente de x**k

    result = S(0)
    for k, c in enumerate(coeffs):
        if c == 0:
            continue
        if k == 0:
            result = result + c
        else:
            result = result + c * (A ** k)
    return result


def qexp(A, order=8, x=None):
    """
    exp(A) truncada a `order` términos de la serie de Taylor.

    IMPORTANTE (rendimiento): si necesitas un valor numérico concreto,
    mantén la constante como un símbolo de sympy dentro de A (ej.
    I*theta*Xgate con theta=var('theta')) y sustituye el valor numérico
    DESPUÉS con .subs(), en vez de pasar un número/complejo ya evaluado
    directamente. Elevar un operador con coeficientes numéricos "sucios"
    (floats/complejos) a potencias altas hace que _evaluar tenga que
    comparar/combinar términos simbólicos sin la limpieza que da trabajar
    con números racionales/simbólicos, y puede volverse muy lento.
    """
    x = x or Symbol('_qseries_x')
    return qseries(sympy_exp(x), A, x, order)


def qsin(A, order=8, x=None):
    """sin(A) truncada a `order` términos de la serie de Taylor."""
    x = x or Symbol('_qseries_x')
    return qseries(sympy_sin(x), A, x, order)


def qcos(A, order=8, x=None):
    """cos(A) truncada a `order` términos de la serie de Taylor."""
    x = x or Symbol('_qseries_x')
    return qseries(sympy_cos(x), A, x, order)