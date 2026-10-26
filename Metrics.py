"""
Port de metrics.py sobre openket_v2.

NOTA IMPORTANTE sobre la diferencia de arquitectura:
En el original, _DiracSum/_DiracMult heredan de DiracObject, así que un
objeto ya "expandido" (_expandir(A)) sigue pudiendo multiplicarse/sumarse
directamente. En openket_v2, _Sum/_Mult son estructuras "scratch" internas
que NO heredan de DiracObject (a propósito, para poder tener resultados
inmutables tipo sympy.Expr). Por eso aquí NUNCA reasigno la variable
pública `A` a su forma scratch cuando todavía la voy a usar en aritmética
(+,-,*) más adelante: uso `_iter_terms_factors(A)` solo para *inspeccionar*
tags/autovalores, dejando `A` intacta para seguir operando con ella.
"""

from sympy import Matrix, I, symbols, Add, sympify

from openket_V2 import (
    Ket, Bra, Operator, AdjointOperator, CreationOperator, AnnihilationOperator,
    DiracMul, DiracAdd, _to_scratch, _finish, _expandir, _Sum, _Mult,
)

__all__ = ['dag', 'comm', 'ptrace', 'trace', 'normalize', 'qmatrix', 'op2dict', 'sub_qexpr']


def _iter_terms_factors(qobj):
    """Devuelve una lista de listas de factores (una lista por término),
    solo para INSPECCIONAR (tags, autovalores). No usar el resultado para
    seguir operando -- para eso, sigue usando el objeto público original."""
    scratch = _expandir(_to_scratch(qobj))
    out = []
    for t in scratch.terms:
        if isinstance(t, _Mult):
            out.append(t.factors)
        else:
            out.append([])  # término puramente escalar
    return out


def dag(qobj):
    """
    Conjugado hermítico de un objeto cuántico (Bra, Ket, Operator, número,
    o una combinación lineal de estos).
    """
    if isinstance(qobj, Bra):
        return Ket(qobj.eig, qobj.op)
    elif isinstance(qobj, Ket):
        return Bra(qobj.eig, qobj.op)
    elif isinstance(qobj, Operator):
        return AdjointOperator(qobj.op)
    elif isinstance(qobj, AdjointOperator):
        # CORREGIDO: el original llamaba a `op(qobj.op)`, una función que no
        # existe en ningún lado del paquete (NameError garantizado). El
        # dagger de un AdjointOperator debe devolver el Operator original.
        return Operator(qobj.op)
    elif isinstance(qobj, DiracAdd):
        return DiracAdd(*[dag(t) for t in qobj.terms])
    elif isinstance(qobj, DiracMul):
        L = list(qobj.factors)
        L.reverse()
        return DiracMul(qobj.coef.conjugate(), *[dag(f) for f in L])
    else:
        return sympify(qobj).conjugate()


def comm(A, B):
    """Conmutador [A, B] = A*B - B*A."""
    return A * B - B * A


def ptrace(A, htag):
    """Traza parcial de A sobre el espacio de Hilbert `htag`."""
    term_factors = _iter_terms_factors(A)
    htags = []
    elements = []
    for factors in term_factors:
        for j in factors:
            htags.append(j.op)
            elements.append(j)
    htags = list(set(htags))
    d = {i: [] for i in htags}
    for i in d:
        for j in elements:
            if j.op == i:
                d[i].append(j.eig)
    for i in d:
        d[i] = list(set(d[i]))
    if htag in d:
        term = 0
        for j in d[htag]:
            term = term + Bra(j, htag) * A * Ket(j, htag)
        A = term
    return A


def trace(A, basis='default'):
    """Traza total de A (todas las Hilbert spaces) o traza en una `basis` dada."""
    if basis == 'default':
        term_factors = _iter_terms_factors(A)
        htags = []
        for factors in term_factors:
            for j in factors:
                htags.append(j.op)
        htags = list(set(htags))
        for i in htags:
            A = ptrace(A, i)
    else:
        temp = 0
        for l in basis:
            temp = temp + dag(l) * A * l
        A = temp
    return A


def normalize(state):
    """Normaliza un vector de estado a norma unidad."""
    norm = dag(state) * state
    norm = float(norm)
    state = state * norm ** -0.5
    return state


def qmatrix(A, basis='default'):
    """Representación matricial de un operador (sympy.Matrix)."""
    scratch = _to_scratch(A)
    scratch = _Sum(scratch) if not isinstance(scratch, _Sum) else scratch
    if basis == 'default':
        eigs = []
        ops = []
        for term in scratch.terms:
            t = term if isinstance(term, _Mult) else _Mult(1, term)
            for fact in t.factors:
                if isinstance(fact, (Ket, Bra)):
                    eigs.append(fact.eig)
                    ops.append(fact.op)
        eigs = sorted(set(eigs), key=str)
        ops = list(set(ops))
        if len(ops) > 1:
            raise Exception('Matrix form can only be obtained for operators '
                             'acting on a single Hilbert space.')
        oper = ops[0]
        M = []
        for i in eigs:
            row = []
            for j in eigs:
                row.append(Bra(i, oper) * A * Ket(j, oper))
            M.append(row)
        M = Matrix(M)
    else:
        M = []
        for i in basis:
            row = []
            for j in basis:
                row.append(dag(i) * A * j)
            M.append(row)
        M = Matrix(M)
    return M


def op2dict(op, basis):
    """Diccionario {<i|op|j> : variable compleja simbólica} para una base dada."""
    n = len(basis)
    re = []
    im = []
    t = []
    for i in range(n * n):
        re.append(symbols('Re%d' % i, real=True))
        im.append(symbols('Im%d' % i, real=True))
    for i in range(len(re)):
        t.append(re[i] + I * im[i])
    D = {}
    k = 0
    for i in range(n):
        for j in range(n):
            D[str(dag(basis[i]) * op * basis[j])] = t[k]
            k += 1
    return D


def sub_qexpr(qexpr, dic):
    """Sustituye en qexpr las claves de `dic` (strings de <i|op|j>) por sus valores."""
    if qexpr == 0.0:
        l = [0]
    else:
        scratch = _to_scratch(qexpr)
        scratch = _Sum(scratch) if not isinstance(scratch, _Sum) else scratch
        l = []
        for i in scratch.terms:
            t = i if isinstance(i, _Mult) else _Mult(1, i)
            y = DiracMul(1, *t.factors)
            key = str(y)
            if key in dic:
                y = dic[key]
            l.append(y * t.coef)
    return Add(*l)