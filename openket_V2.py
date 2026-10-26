"""
OpenKet v2 -- port fiel del motor original (_evaluar/_expandir/_canonic_mult/
_canonic_suma/_is_equal) integrado a sympy.

DISEÑO:
- Las hojas (Ket, Bra, Operator, AdjointOperator, CreationOperator,
  AnnihilationOperator) son sympy.Expr inmutables: se pueden usar en subs(),
  imprimir con el printer de sympy, comparar estructuralmente, etc.
- El motor de evaluación (_evaluar y compañía) NO cambia de algoritmo: sigue
  operando sobre listas mutables (_Sum/_Mult), exactamente como el original,
  porque esa es la única forma sensata de portar un algoritmo de escaneo
  hacia adelante/atrás a través de una cadena de factores. La diferencia es
  que estas son ahora clases "scratch" internas (no Expr, no se exponen).
- Al final de cada operación (+,-,*,/,**) el resultado scratch se convierte
  ("_finish") a un objeto inmutable público: DiracMul o DiracAdd (o un
  escalar plano de sympy si todo se cancela).
- _op_priority alto es lo que permite que expresiones como sqrt(n)*Ket(m)
  (donde el operando izquierdo es un Expr genérico de sympy) le cedan el
  paso a nuestra lógica en vez de construir un Mul genérico de sympy.
"""

from sympy import Expr, Add, sympify, S, Integer
from sympy import integrate as _sympy_integrate

__all__ = ["Ket", "Bra", "Operator", "AdjointOperator",
           "CreationOperator", "AnnihilationOperator",
           "DiracMul", "DiracAdd", "qintegrate",
           "_to_scratch", "_finish", "_expandir", "_Sum", "_Mult"]


# =====================================================================
# Hojas públicas: Ket, Bra, Operator, AdjointOperator, Creation/Annihilation
# =====================================================================

class DiracObject(Expr):
    is_commutative = False
    _op_priority = 20.0   # más alto que el default de Expr (~10-11):
                           # obliga a que sympy nos ceda el paso incluso
                           # cuando NO somos el operando izquierdo.

    def __add__(self, other):
        return _finish(_evaluar(_Sum(_to_scratch(self), _to_scratch(other))))

    def __radd__(self, other):
        return _finish(_evaluar(_Sum(_to_scratch(other), _to_scratch(self))))

    def __sub__(self, other):
        return _finish(_evaluar(_Sum(_to_scratch(self),
                                      _Mult(-1, _to_scratch(other)))))

    def __rsub__(self, other):
        return _finish(_evaluar(_Sum(_to_scratch(other),
                                      _Mult(-1, _to_scratch(self)))))

    def __mul__(self, other):
        a = _to_scratch(self)
        b = _to_scratch(other)
        if isinstance(a, _Mult):
            if isinstance(b, _Mult):
                return _finish(_evaluar(_Mult(a.coef * b.coef,
                                               *(a.factors + b.factors))))
            return _finish(_evaluar(_Mult(a.coef, *(a.factors + [b]))))
        else:
            if isinstance(b, _Mult):
                return _finish(_evaluar(_Mult(b.coef, *([a] + b.factors))))
            return _finish(_evaluar(_Mult(1, a, b)))

    def __rmul__(self, other):
        a = _to_scratch(self)
        b = _to_scratch(other)
        if isinstance(a, _Mult):
            if isinstance(b, _Mult):
                return _finish(_evaluar(_Mult(b.coef * a.coef,
                                               *(b.factors + a.factors))))
            return _finish(_evaluar(_Mult(a.coef, *([b] + a.factors))))
        else:
            if isinstance(b, _Mult):
                return _finish(_evaluar(_Mult(b.coef, *(b.factors + [a]))))
            return _finish(_evaluar(_Mult(1, b, a)))

    def __truediv__(self, other):
        a = _to_scratch(self)
        if isinstance(a, _Mult):
            return _finish(_evaluar(_Mult(a.coef / other, *a.factors)))
        return _finish(_evaluar(_Mult(Add(1) / other, a)))

    def __pow__(self, n):
        # IMPORTANTE (rendimiento): multiplicar de forma incremental, NO
        # construir _Mult(1, self, self, ..., self) de una sola vez. Si self
        # es una suma (_Sum) de k términos, meterla n veces en un solo _Mult
        # obliga a _expandir a distribuir TODO antes de poder colapsar nada
        # por ortogonalidad -- eso es una explosión combinatoria de k**n
        # términos (verificado: cuadruplicaba el tiempo en cada orden extra
        # para una suma de 2 términos). Multiplicando paso a paso, cada
        # producto se colapsa de inmediato antes de seguir, evitando la
        # explosión.
        if n == 0:
            return S(1)
        result = self
        for _ in range(n - 1):
            result = result * self
        return result


class Ket(DiracObject):
    """Ket(a, A). API idéntica al original: .eig, .op, repr '|a>' / '|a_A>'."""
    def __new__(cls, a, A='default'):
        eig = sympify(a)
        obj = Expr.__new__(cls, eig, S(A) if isinstance(A, str) else sympify(A))
        obj.eig = eig
        obj.op = A
        return obj

    def __repr__(self):
        if self.op == 'default':
            return "|" + repr(self.eig) + ">"
        return "|" + str(self.eig) + "_" + str(self.op) + ">"

    def _sympystr(self, printer):
        return self.__repr__()


class Bra(DiracObject):
    """Bra(a, A). API idéntica al original: .eig, .op, repr '<a|' / '<a_A|'."""
    def __new__(cls, a, A='default'):
        eig = sympify(a)
        obj = Expr.__new__(cls, eig, S(A) if isinstance(A, str) else sympify(A))
        obj.eig = eig
        obj.op = A
        return obj

    def __repr__(self):
        if self.op == 'default':
            return "<" + str(self.eig) + "|"
        return "<" + str(self.eig) + "_" + str(self.op) + "|"

    def _sympystr(self, printer):
        return self.__repr__()


class Operator(DiracObject):
    """Operator(A). API idéntica al original: .op, repr 'A'."""
    def __new__(cls, A):
        obj = Expr.__new__(cls, S(A) if isinstance(A, str) else sympify(A))
        obj.op = A
        return obj

    def __repr__(self):
        return str(self.op)

    def _sympystr(self, printer):
        return str(self.op)


class AdjointOperator(DiracObject):
    """Dual del operador. API idéntica al original."""
    def __new__(cls, A):
        obj = Expr.__new__(cls, S(A) if isinstance(A, str) else sympify(A))
        obj.op = A
        return obj

    def __repr__(self):
        return str(self.op) + "'"

    def _sympystr(self, printer):
        return self.__repr__()


class CreationOperator(DiracObject):
    r"""
    a^+|n> = sqrt(n+1)|n+1>      <n|a^+ = sqrt(n)<n-1|
    htag='default' por defecto, maximum=0 significa "sin cota".
    """
    def __new__(cls, htag='default', maximum=0):
        obj = Expr.__new__(cls, S(htag) if isinstance(htag, str) else sympify(htag),
                            sympify(maximum))
        obj.op = htag
        obj.maximum = maximum
        return obj

    def __repr__(self):
        return "a'" if self.op == 'default' else "a'(" + str(self.op) + ")"

    def _sympystr(self, printer):
        return self.__repr__()


class AnnihilationOperator(DiracObject):
    r"""
    a|n> = sqrt(n)|n-1>          <n|a = sqrt(n+1)<n+1|
    htag='default' por defecto, maximum=0 significa "sin cota".
    """
    def __new__(cls, htag='default', maximum=0):
        obj = Expr.__new__(cls, S(htag) if isinstance(htag, str) else sympify(htag),
                            sympify(maximum))
        obj.op = htag
        obj.maximum = maximum
        return obj

    def __repr__(self):
        return "a" if self.op == 'default' else "a(" + str(self.op) + ")"

    def _sympystr(self, printer):
        return self.__repr__()


# =====================================================================
# Resultados públicos inmutables (lo que ve el usuario tras +,-,*,/,**)
# =====================================================================

class DiracMul(DiracObject):
    """
    Resultado inmutable de una 'multiplicación' evaluada: un coeficiente
    (self.coef, sympy) y una lista de factores restantes sin colapsar
    (self.factors: Ket/Bra/Operator/... que no se pudieron simplificar más).
    Equivalente inmutable del _DiracMult original.
    """
    def __new__(cls, coef, *factors):
        coef = sympify(coef)
        obj = Expr.__new__(cls, coef, *factors)
        obj.coef = coef
        obj.factors = list(factors)
        return obj

    def __repr__(self):
        if self.coef.is_Add:
            s = '(' + str(self.coef) + ')'
        elif self.coef == 1:
            s = ''
        elif self.coef == -1:
            s = ' - '
        else:
            atomos = list(self.coef.atoms())
            if atomos and atomos[0].is_Number:
                s = ' - ' + str(-self.coef) if atomos[0] < 0 else str(self.coef)
            else:
                s = str(self.coef)
        n = 0
        k = 1
        while n <= len(self.factors) - 1:
            if n == len(self.factors) - 1:
                s += str(self.factors[n]) + '**' + str(k) if k > 1 else str(self.factors[n])
                break
            if _is_equal(self.factors[n], self.factors[n + 1]):
                k += 1
                n += 1
                continue
            if k > 1:
                s += '(' + str(self.factors[n]) + '**' + str(k) + ')'
                k = 1
            else:
                s += str(self.factors[n])
            n += 1
        return s

    def _sympystr(self, printer):
        return self.__repr__()

    def _eval_derivative(self, x):
        """Regla del producto, reutilizando nuestros propios operadores +,*
        (que ya sabemos que simplifican/combinan correctamente). Cada factor
        se diferencia con su propio .diff(x): si no depende de x, sympy ya
        devuelve 0 automáticamente (heredado de Expr); si depende de x pero
        no sabemos cómo diferenciarlo físicamente (p.ej. un Ket(x) genérico),
        queda como Derivative(...) sin evaluar -- honesto, no inventamos una
        regla física que no existe de forma genérica."""
        parts = list(self.args)  # (coef, factor1, factor2, ...)
        result = S(0)
        for i in range(len(parts)):
            d = parts[i].diff(x)
            if d == 0:
                continue
            new_parts = list(parts)
            new_parts[i] = d
            term = new_parts[0]
            for p in new_parts[1:]:
                term = term * p
            result = result + term
        return result


class DiracAdd(DiracObject):
    """Suma inmutable de términos (escalares y/o DiracMul). Equivalente
    inmutable del _DiracSum original."""
    def __new__(cls, *terms):
        obj = Expr.__new__(cls, *terms)
        obj.terms = list(terms)
        return obj

    def __repr__(self):
        s = str(self.terms[0])
        for term in self.terms[1:]:
            term_str = str(term)
            if term_str.startswith(' - '):
                s += term_str
            else:
                s += ' + ' + term_str
        return s

    def _sympystr(self, printer):
        return self.__repr__()

    def _eval_derivative(self, x):
        """Regla de la suma: derivada de cada término, sumadas con nuestro
        propio '+' (así se simplifican/combinan igual que cualquier otra
        suma de nuestra álgebra)."""
        result = S(0)
        for t in self.terms:
            result = result + t.diff(x)
        return result


# =====================================================================
# Puente scratch <-> público
# =====================================================================

def _to_scratch(x):
    """Convierte lo que sea (leaf, DiracMul, DiracAdd, escalar) a la
    representación mutable interna que usa el algoritmo original."""
    if isinstance(x, DiracAdd):
        return _Sum(*[_to_scratch(t) for t in x.terms])
    if isinstance(x, DiracMul):
        return _Mult(x.coef, *x.factors)
    return x  # leaf (Ket/Bra/Operator/...) o escalar: se usan tal cual,
              # igual que en el original.


def _finish(x):
    """Convierte el resultado final de _evaluar (scratch) al tipo público
    inmutable correspondiente."""
    if isinstance(x, _Sum):
        return DiracAdd(*[_finish(t) for t in x.terms])
    if isinstance(x, _Mult):
        return DiracMul(x.coef, *x.factors)
    return sympify(x)


def qintegrate(expr, *args, **kwargs):
    """
    Integra una expresión de OpenKet respecto a una variable ordinaria
    (no un espacio de Hilbert). A diferencia de diff(), sympy.integrate()
    NO tiene un gancho generico (_eval_derivative) para objetos custom --
    sus algoritmos de integracion reconocen literalmente sympy.Mul/Add,
    no DiracMul/DiracAdd. Por eso esta funcion existe aparte.

    LIMITACION IMPORTANTE (explicita, no oculta): esto solo funciona
    correctamente cuando los Ket/Bra/Operator involucrados NO dependen de
    la variable de integracion (el caso comun: f(x)*Ket(n) con n fijo).
    Si algun factor SI depende de la variable de integracion (ej. Ket(x)),
    se lanza un error en vez de dar un resultado silenciosamente incorrecto
    -- integrar sobre un ket que varia continuamente (resolucion de la
    identidad, representacion de posicion, etc.) requiere reglas fisicas
    especificas que no se pueden asumir de forma generica.
    """
    var = args[0] if args and not isinstance(args[0], (tuple, list)) else (args[0][0] if args else None)

    def _integrate_leaf(coef, factors):
        for f in factors:
            if var is not None and var in f.free_symbols:
                raise NotImplementedError(
                    f"qintegrate: el factor {f} depende de la variable de "
                    f"integracion {var}. Esto requiere una regla fisica "
                    f"explicita (p.ej. representacion de posicion) que "
                    f"qintegrate no puede asumir de forma generica."
                )
        new_coef = _sympy_integrate(coef, *args, **kwargs)
        if len(factors) == 0:
            return new_coef
        result = new_coef
        for f in factors:
            result = result * f
        return result

    if isinstance(expr, DiracAdd):
        result = S(0)
        for t in expr.terms:
            result = result + qintegrate(t, *args, **kwargs)
        return result
    if isinstance(expr, DiracMul):
        return _integrate_leaf(expr.coef, expr.factors)
    if isinstance(expr, (Ket, Bra, Operator, AdjointOperator, CreationOperator, AnnihilationOperator)):
        return _integrate_leaf(S(1), [expr])
    # expresión sin nada de OpenKet: delega directo a sympy
    return _sympy_integrate(expr, *args, **kwargs)


# =====================================================================
# Motor interno (scratch, mutable) -- port ~1:1 del original
# =====================================================================

class _Sum(object):
    """Equivalente mutable/interno de _DiracSum."""
    def __init__(self, *terms):
        self.terms = list(terms)
        for term in self.terms:
            if isinstance(term, _Sum):
                n = self.terms.index(term)
                self.terms[n:n + 1] = term.terms


class _Mult(object):
    """Equivalente mutable/interno de _DiracMult."""
    def __init__(self, coefficient, *factors):
        self.coef = Add(coefficient)
        self.factors = list(factors)
        for factor in self.factors:
            if isinstance(factor, _Mult):
                self.coef = self.coef * factor.coef
                n = self.factors.index(factor)
                self.factors[n:n + 1] = factor.factors


def _evaluar(expression):
    """Idéntico en algoritmo al _evaluar original: expande, aplica 
    operadores buscando hacia adelante/atrás en la cadena de factores, 
    colapsa productos internos, ordena en forma canónica y 
    combina términos iguales."""
    expression = _expandir(expression)

    for n in range(len(expression.terms)):
        if not isinstance(expression.terms[n], _Mult):
            expression.terms[n] = _Mult(1, expression.terms[n])

    for term in expression.terms:
        for n in range(len(term.factors)):
            if isinstance(term.factors[n], AdjointOperator):
                for m in reversed(range(n)):
                    if isinstance(term.factors[m], Bra):
                        if term.factors[n].op == term.factors[m].op:
                            term.coef = term.coef * term.factors[m].eig.conjugate()
                            term.factors[n] = 'out'
                            break
                    elif isinstance(term.factors[m], (Operator, AdjointOperator,
                                                        CreationOperator, AnnihilationOperator)):
                        break

        for n in reversed(range(len(term.factors))):
            if isinstance(term.factors[n], Ket):
                pass
            elif isinstance(term.factors[n], Bra):
                for m in range(n + 1, len(term.factors)):
                    if isinstance(term.factors[m], Ket):
                        if term.factors[n].op == term.factors[m].op:
                            if term.factors[n].eig == term.factors[m].eig:
                                term.factors[n] = 'out'
                                term.factors[m] = 'out'
                            else:
                                term.coef = 0
                            break
                    elif isinstance(term.factors[m], (Operator, AdjointOperator)):
                        break
            elif isinstance(term.factors[n], Operator):
                for m in range(n + 1, len(term.factors)):
                    if isinstance(term.factors[m], Ket):
                        if term.factors[n].op == term.factors[m].op:
                            term.coef = term.coef * term.factors[m].eig
                            term.factors[n] = 'out'
                            break
                    elif isinstance(term.factors[m], (Operator, AdjointOperator,
                                                        CreationOperator, AnnihilationOperator)):
                        break
            elif isinstance(term.factors[n], AdjointOperator):
                pass
            elif isinstance(term.factors[n], CreationOperator):
                for m in range(n + 1, len(term.factors)):
                    if isinstance(term.factors[m], Ket):
                        if term.factors[n].op == term.factors[m].op:
                            if term.factors[n].maximum == 0:
                                term.coef = term.coef * (term.factors[m].eig + 1) ** (Add(1) / Add(2))
                                term.factors[m] = Ket(term.factors[m].eig + 1, term.factors[m].op)
                                term.factors[n] = 'out'
                            else:
                                if term.factors[m].eig < term.factors[n].maximum:
                                    term.coef = term.coef * (term.factors[m].eig + 1) ** (Add(1) / Add(2))
                                    term.factors[m] = Ket(term.factors[m].eig + 1, term.factors[m].op)
                                    term.factors[n] = 'out'
                                else:
                                    term.coef = Add(0)
                            break
                    elif isinstance(term.factors[m], (Operator, AdjointOperator,
                                                        CreationOperator, AnnihilationOperator)):
                        break
                if term.factors[n] != 'out':
                    for m in reversed(range(n)):
                        if isinstance(term.factors[m], Bra):
                            if term.factors[n].op == term.factors[m].op:
                                term.coef = term.coef * (term.factors[m].eig) ** (Add(1) / Add(2))
                                term.factors[m] = Bra(term.factors[m].eig - 1, term.factors[m].op)
                                term.factors[n] = 'out'
                                break
                        elif isinstance(term.factors[m], (Operator, AdjointOperator,
                                                            CreationOperator, AnnihilationOperator)):
                            break
            elif isinstance(term.factors[n], AnnihilationOperator):
                for m in range(n + 1, len(term.factors)):
                    if isinstance(term.factors[m], Ket):
                        if term.factors[n].op == term.factors[m].op:
                            term.coef = term.coef * (term.factors[m].eig) ** (Add(1) / Add(2))
                            term.factors[m] = Ket(term.factors[m].eig - 1, term.factors[m].op)
                            term.factors[n] = 'out'
                            break
                    elif isinstance(term.factors[m], (Operator, AdjointOperator,
                                                        CreationOperator, AnnihilationOperator)):
                        break
                if term.factors[n] != 'out':
                    for m in reversed(range(n)):
                        if isinstance(term.factors[m], Bra):
                            if term.factors[n].op == term.factors[m].op:
                                if term.factors[n].maximum == 0:
                                    term.coef = term.coef * (term.factors[m].eig + 1) ** (Add(1) / Add(2))
                                    term.factors[m] = Bra(term.factors[m].eig + 1, term.factors[m].op)
                                    term.factors[n] = 'out'
                                else:
                                    if term.factors[m].eig < term.factors[n].maximum:
                                        term.coef = term.coef * (term.factors[m].eig + 1) ** (Add(1) / Add(2))
                                        term.factors[m] = Bra(term.factors[m].eig + 1, term.factors[m].op)
                                        term.factors[n] = 'out'
                                    else:
                                        term.coef = Add(0)
                                break
                        elif isinstance(term.factors[m], (Operator, AdjointOperator,
                                                            CreationOperator, AnnihilationOperator)):
                            break
            else:
                if term.factors[n] != 'out':
                    term.coef = term.coef * term.factors[n]
                    term.factors[n] = 'out'
        while 'out' in term.factors:
            term.factors.remove('out')

    for n in range(len(expression.terms)):
        expression.terms[n] = _canonic_mult(expression.terms[n])

    for n in range(len(expression.terms)):
        for m in range(n + 1, len(expression.terms)):
            if _is_equal(expression.terms[n], expression.terms[m]):
                expression.terms[n] = _Mult(expression.terms[n].coef + expression.terms[m].coef,
                                             *expression.terms[n].factors)
                expression.terms[m] = 'out'
    while 'out' in expression.terms:
        expression.terms.remove('out')

    for n in range(len(expression.terms)):
        if expression.terms[n].coef == 0:
            expression.terms[n] = 'out'
    while 'out' in expression.terms:
        expression.terms.remove('out')

    for n in range(len(expression.terms)):
        if len(expression.terms[n].factors) == 0:
            expression.terms[n] = expression.terms[n].coef

    expression = _canonic_suma(expression)
    return expression


def _expandir(expression):
    if isinstance(expression, _Sum):
        nested = True
        while nested:
            n = 0
            nmax = len(expression.terms)
            while n <= nmax:
                if n == nmax:
                    nested = False
                    break
                if isinstance(expression.terms[n], _Mult):
                    k = 0
                    kmax = len(expression.terms[n].factors)
                    while k < kmax:
                        if isinstance(expression.terms[n].factors[k], _Sum):
                            L = []
                            for term in expression.terms[n].factors[k].terms:
                                L.append(_Mult(expression.terms[n].coef,
                                                *expression.terms[n].factors[:k] + [term]
                                                + expression.terms[n].factors[k + 1:]))
                            expression = _Sum(*expression.terms[:n] + L + expression.terms[n + 1:])
                            k = kmax
                            n = nmax
                        else:
                            k += 1
                    n += 1
                else:
                    n += 1
        return expression
    elif isinstance(expression, _Mult):
        return _expandir(_Sum(expression))
    else:
        return _Sum(expression)


def _canonic_mult(term):
    L = []
    kets = []
    bras = []
    position_lock = True
    for k, u in enumerate(term.factors):
        if isinstance(u, (Operator, AdjointOperator, CreationOperator, AnnihilationOperator)):
            if not position_lock:
                kets.sort()
                bras.sort()
                bras.reverse()
                for t in kets:
                    L.append(term.factors[t[2]])
                for t in bras:
                    L.append(term.factors[t[2]])
                L.append(u)
                kets = []
                bras = []
                position_lock = True
            else:
                L.append(u)
        elif isinstance(u, Ket):
            position_lock = False
            kets.append((str(u.op), u.eig, k))
        elif isinstance(u, Bra):
            position_lock = False
            bras.append((str(u.op), u.eig, k))
        else:
            print('Error in canonic_mult. Type not valid.')
    if len(kets) != 0 or len(bras) != 0:
        kets.sort()
        bras.sort()
        bras.reverse()
        for t in kets:
            L.append(term.factors[t[2]])
        for t in bras:
            L.append(term.factors[t[2]])
    term = _Mult(term.coef, *L)
    return term


def _canonic_suma(expression):
    L1 = []
    indep = 0
    for k, u in enumerate(expression.terms):
        if isinstance(u, _Mult):
            if isinstance(u.factors[0], CreationOperator):
                L1.append((1, str(u.factors[0].op), 0, k))
            if isinstance(u.factors[0], AnnihilationOperator):
                L1.append((2, str(u.factors[0].op), 0, k))
            if isinstance(u.factors[0], AdjointOperator):
                L1.append((3, str(u.factors[0].op), 0, k))
            elif isinstance(u.factors[0], Operator):
                L1.append((4, str(u.factors[0].op), 0, k))
            elif isinstance(u.factors[0], Ket):
                L1.append((5, str(u.factors[0].op), u.factors[0].eig, k))
            elif isinstance(u.factors[0], Bra):
                L1.append((6, str(u.factors[0].op), u.factors[0].eig, k))
        else:
            indep = indep + u
    L1.sort(key=lambda t: (t[0], t[1], str(t[2])))
    L2 = []
    if indep != 0:
        L2.append(indep)
    for t in L1:
        L2.append(expression.terms[t[3]])
    if len(L2) == 0:
        expression = 0
    elif len(L2) == 1:
        expression = L2[0]
    else:
        expression = _Sum(*L2)
    return expression


def _is_equal(u, v):
    if type(u) != type(v):
        return False
    if u == 'out':
        return False
    if isinstance(u, (Ket, Bra)):
        return u.op == v.op and u.eig == v.eig
    if isinstance(u, (Operator, AdjointOperator, CreationOperator, AnnihilationOperator)):
        return u.op == v.op
    if isinstance(u, _Mult):
        if len(u.factors) != len(v.factors):
            return False
        for k in range(len(u.factors)):
            if not _is_equal(u.factors[k], v.factors[k]):
                return False
        return True
    return False