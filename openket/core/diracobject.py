from sympy import Add

class DiracObject(object):
    """
    A class for representing quantum objects, such as quantum operators and states.

    The DiracObject class is the OpenKet representation of quantum operators and state vectors.
    This class also implements math operations +,-,* between DiracObject instances (and / by a complex number).
    """
    def __init__(self):
        pass
    def __add__(self, other):
        return _evaluar(DiracSum(self, other))
    def __radd__(self, other):
        return _evaluar(DiracSum(other, self))
    def __sub__(self, other):
        return _evaluar(DiracSum(self, DiracMult(-1, other)))
    def __rsub__(self, other):
        return _evaluar(DiracSum(other, DiracMult(-1, self)))
    def __mul__(self, other):
        if isinstance(self, DiracMult):
            if isinstance(other, DiracMult):
                return _evaluar(DiracMult(self.coef*other.coef, *self.factors \
                                    + other.factors))
            else:
                return _evaluar(DiracMult(self.coef, *self.factors + [other]))
        else:
            if isinstance(other, DiracMult):
                return _evaluar(DiracMult(other.coef, *[self] + other.factors))
            else:
                return _evaluar(DiracMult(1, self, other))
    def __rmul__(self, other):
        if isinstance(self, DiracMult):
            if isinstance(other, DiracMult):
                return _evaluar(DiracMult(other.coef*self.coef, *other.factors \
                                    + self.factors))
            else:
                return _evaluar(DiracMult(self.coef, *[other] + self.factors))
        else:
            if isinstance(other, DiracMult):
                return _evaluar(DiracMult(other.coef, *other.factors + [self]))
            else:
                return _evaluar(DiracMult(1, other, self))
    def __div__(self, other):
        if isinstance(self, DiracMult):
            return _evaluar(DiracMult(self.coef/other, *self.factors))
        else:
            return _evaluar(DiracMult(Add(1)/other, self))
    def __pow__(self, n):
        L = []
        for k in range(n):
            L.append(self)
        return _evaluar(DiracMult(1, *L))
    
class Ket(DiracObject):
    """
    A class for representing a state, called ket. `Ket` is a class inheritanced from `DiracObject`.

    :param a: Vector tag. Must be either a number class or a Sympy variable.
    :type a: int or Sympy variable
    :param A: Operator tag. The ket is then the eigenvector of the operator, and they comply with the eigenvalues equation.
                Also this param could be the Hilbert space tag where the ket lives.
    :type A: string, optional
    
    """
    def __init__(self, a, A = 'default'):
        self.eig = Add(a)   # This is just a way to make any number compatible
        self.op = A         # with Sympy functions.
    def __repr__(self):
        if self.op == 'default':
            return "|" + repr(self.eig) + ">"
        else:
            return "|" + str(self.eig) + "_" + str(self.op) + ">"

class Bra(DiracObject):
    """
    A class for representing the vectors corresponding to the dual space, called bras.
    `Bra` is a class inheritanced from `DiracObject`.

    :param a: Vector dual tag. Must be either a number class or a Sympy variable.
    :type a: int or Sympy variable
    :param A: Operator tag representing the eigenvalues equation for vector whose bra is the dual.
                Also this param could be the Hilbert space tag where the corresponding ket lives.
    :type A: string, optional

    """
    def __init__(self, a, A = 'default'):
        self.eig = Add(a)
        self.op = A
    def __repr__(self):
        if self.op == 'default':
            return "<" + str(self.eig) + "|"
        else:
            return "<" + str(self.eig) + "_" + str(self.op) + "|"

class Operator(DiracObject):
    """
    A class for representing an operator.
    `Operator` is a class inheritanced from `DiracObject`.

    :param A: Operator tag.
    :type A: string

    """
    def __init__(self, A):
        self.op = A
    def __repr__(self):
        return str(self.op)

class AdjointOperator(DiracObject):
    """Definition of the dual operator class."""
    def __init__(self, A):
        self.op = A
    def __repr__(self):
        return str(self.op) + "'"

class CreationOperator(DiracObject):
    """
    A class for representing the Creation operator.

    .. math::

        \\begin{aligned}
        \hat{a}^{\dagger}|n\\rangle &= \sqrt{n+1}|n+1\\rangle \\\\
        \\langle n| \hat{a}^{\dagger} &= \sqrt{n} \\langle n-1| \\\\
        \\end{aligned}

    :param htag: Hilbert space tag.
    :type htag: string, optional
    :param maximum: Maximum value for n.
    :type maximum: int, optional

    """
    def __init__(self, htag = 'default', maximum = 0):
        self.op = htag
        self.maximum = maximum
    def __repr__(self):
        if self.op == 'default':
            return "a'"
        else:
            return "a'(" + str(self.op) + ")"

class AnnihilationOperator(DiracObject):
    """
    A class for representing the Annihilation operator.

    .. math::

        \\begin{aligned}
        \hat{a}|n\\rangle &= \sqrt{n}|n-1\\rangle \\\\
        \\langle n| \hat{a} &= \sqrt{n+1} \\langle n+1| \\\\
        \\end{aligned}

    :param htag: Hilbert space tag.
    :type htag: string, optional
    :param maximum: Maximum value for n.
    :type maximum: int, optional

    """
    def __init__(self, htag = 'default', maximum = 0):
        self.op = htag
        self.maximum = maximum
    def __repr__(self):
        if self.op == 'default':
            return 'a'
        else:
            return 'a(' + str(self.op) + ')'

class DiracSum(DiracObject):
    """This class represents the sum of several elements, which
    can be kets, bras, operators, dual operators, numbers, sums
    and 'multiplications', and stores them in a list."""
    def __init__(self, *terms):
        self.terms = list(terms)
        for term in self.terms:
            if isinstance(term, DiracSum):
                n = self.terms.index(term)
                self.terms[n:n+1] = term.terms
    def __repr__(self):
        s = str(self.terms[0])
        for term in self.terms[1:]:
            # if isinstance(term, DiracMult):                    #It is not longer needed
            #     if sorted(list(term.coef.atoms()))[0] < 0:     #to correctly display the sum as sympy changed 
            #         s += str(term)
            #     else:
            #         s += ' + ' + str(term)
            # else:
                s += ' + ' + str(term)
        return s

class DiracMult(DiracObject):
    """This class represents the 'multiplication' of any number of
    elements, which can be kets, bras, operators, dual operators,
    scalars, sums (DiracSum) and even other 'multiplications' (DiracMult),
    and stores them in a list. This 'multiplication' can have different
    meanings depending on the arguments: inner product, outer product,
    tensor product, application of operators, or plain multiplication.
    It also has an attribute (self.coef) which contains
    the multiplicative coefficient. When called
    directly (i.e., not by overloading the operator *), a multiplicative
    coefficient must always be included as the first argument. If there
    is no such coefficient, 1 must be the first argument."""
    def __init__(self, coefficient, *factors):
        self.coef = Add(coefficient)
        self.factors = list(factors)
        for factor in self.factors:
            if isinstance(factor, DiracMult):
                self.coef = self.coef*factor.coef
                n = self.factors.index(factor)
                self.factors[n:n+1] = factor.factors
    def __repr__(self):
        if self.coef.is_Add:
            s = '(' + str(self.coef) + ')'
        elif self.coef == 1:
            s = ''
        elif self.coef == -1:
            s = ' - '
        else:
            atomos = list(self.coef.atoms()) # Correct: It would be better to not transform to list
            if atomos[0].is_Number:
                if atomos[0] < 0:
                    s = ' - ' + str(-self.coef)
                else:
                    s = str(self.coef)
            else:
                s = str(self.coef)
        """The following lines allow expressions such as
        AAAB to be printed as (A**3)B."""
        n = 0; k = 1    # k counts how many times an element is repeated.
        while n <= len(self.factors) - 1:
            if n == len(self.factors) - 1:
                if k > 1:
                    s = s + str(self.factors[n]) + '**' + str(k)
                else:
                    s = s + str(self.factors[n])
                break
            else:
                if _is_equal(self.factors[n], self.factors[n+1]):
                    k = k + 1
                    n = n + 1
                    continue
                else:
                    if k > 1:
                        s = s + '(' + str(self.factors[n]) + '**' + str(k) + ')'
                        k = 1
                        n = n + 1
                    else:
                        s = s + str(self.factors[n])
                        n = n + 1
        return s


def _evaluar(expression):
    """This function evaluates any expression, be it a sum (class DiracSum),
    a 'multiplication' (class DiracMult), or simple objects such as kets (in this
    case the function returns the object without changes). It basically
    works the following way: it first expands the expression by means of the
    function expandir, which always returns a sum object (DiracSum), even
    if it has only one element. The function then operates on each term,
    performing inner products and applying operators
    when possible. Afterwards, each term is sorted into the
    'canonic form' by means of the function _canonic_mult. The next stage
    is to check which terms are 'equal', so they can be added or
    substracted. This is accomplished with help of the function
    is_equal. Finally, the sum terms are sorted into the 'canonic form'
    as defined by the function _canonic_suma."""
    expression = _expandir(expression)
    """The expression is now a single sum of terms, none of which can be
    expanded further. Next, each term is turned into a 'multiplication'."""
    for n in range(len(expression.terms)):
        if isinstance(expression.terms[n], DiracMult):
            pass
        else:
            expression.terms[n] = DiracMult(1, expression.terms[n])
    # Now each term is evaluated:
    for term in expression.terms:
        for n in range(len(term.factors)):
            if isinstance(term.factors[n], AdjointOperator):
                """A dual operator 'looks backwards' for its eigenvector
                (a Bra in this case). If it finds another operator, it
                stops."""
                for m in reversed(range(n)):
                    if isinstance(term.factors[m], Bra):
                        if term.factors[n].op == term.factors[m].op:
                            term.coef = term.coef*term.factors[m].eig.conjugate()
                            term.factors[n] = 'out'
                            break
                    elif isinstance(term.factors[m], Operator) or \
                  isinstance(term.factors[m], AdjointOperator) or isinstance(term.factors[m], CreationOperator) or \
                  isinstance(term.factors[m], AnnihilationOperator):
                        break
        for n in reversed(range(len(term.factors))):
            if isinstance(term.factors[n], Ket):
                pass
            elif isinstance(term.factors[n], Bra):
                for m in range(n+1, len(term.factors)):
                    """Each Bra 'looks forward' for a Ket with the same
                    operator to perform an inner product. If it encounters
                    an operator, it stops."""
                    if isinstance(term.factors[m], Ket):
                        if term.factors[n].op == term.factors[m].op:
                            if term.factors[n].eig == term.factors[m].eig:
                                """The preferred way of eliminating elements
                                without affecting the length of the list is
                                labeling them with the string 'out' and
                                removing them later."""
                                term.factors[n] = 'out'
                                term.factors[m] = 'out'
                            else:
                                term.coef = 0
                            break
                    elif isinstance(term.factors[m], Operator) or \
                           isinstance(term.factors[m], AdjointOperator):
                        break
            elif isinstance(term.factors[n], Operator):
                """An operator 'looks forward' for its eigenvector. If it
                finds another operator, it stops."""
                for m in range(n+1, len(term.factors)):
                    if isinstance(term.factors[m], Ket):
                        if term.factors[n].op == term.factors[m].op:
                            term.coef = term.coef*term.factors[m].eig
                            term.factors[n] = 'out'
                            break
                    elif isinstance(term.factors[m], Operator) or \
                  isinstance(term.factors[m], AdjointOperator) or isinstance(term.factors[m], CreationOperator) or \
                  isinstance(term.factors[m], AnnihilationOperator):
                        break                     
            elif isinstance(term.factors[n], AdjointOperator):
                pass
            elif isinstance(term.factors[n], CreationOperator):
                """A creation operator 'searches forward' for a Ket living in
                the same Hilbert space (with the same operator, in this case)
                to perform a*|n> = sqrt(n+1)|n+1> and 'searches backwards' to
                perform <n|a* = sqrt(n)<n-1|. If it finds another operator,
                it stops. Furthermore, if |n+1> exceeds a maximum value
                (for example, the number of possible states in a system),
                the result is zero."""
                for m in range(n+1, len(term.factors)):
                    if isinstance(term.factors[m], Ket):
                        if term.factors[n].op == term.factors[m].op:
                            if term.factors[n].maximum == 0:
                                term.coef = term.coef*(term.factors[m].eig + 1)**(Add(1)/Add(2))
                                term.factors[m] = Ket(term.factors[m].eig + 1, term.factors[m].op)
                                term.factors[n] = 'out'
                            else:
                                if term.factors[m].eig < term.factors[n].maximum:
                                    term.coef = term.coef*(term.factors[m].eig + 1)**(Add(1)/Add(2))
                                    term.factors[m] = Ket(term.factors[m].eig + 1, term.factors[m].op)
                                    term.factors[n] = 'out'
                                else:
                                    term.coef = Add(0)
                            break
                    elif isinstance(term.factors[m], Operator) or \
                  isinstance(term.factors[m], AdjointOperator) or isinstance(term.factors[m], CreationOperator) or \
                  isinstance(term.factors[m], AnnihilationOperator):
                        break
                if term.factors[n] != 'out':
                    for m in reversed(range(n)):
                        if isinstance(term.factors[m], Bra):
                            if term.factors[n].op == term.factors[m].op:
                                term.coef = term.coef*(term.factors[m].eig)**(Add(1)/Add(2))
                                term.factors[m] = Bra(term.factors[m].eig - 1, term.factors[m].op)
                                term.factors[n] = 'out'
                                break
                        elif isinstance(term.factors[m], Operator) or \
                    isinstance(term.factors[m], AdjointOperator) or isinstance(term.factors[m], CreationOperator) or \
                    isinstance(term.factors[m], AnnihilationOperator):
                            break
            elif isinstance(term.factors[n], AnnihilationOperator):
                """An annihilation operator 'searches forward' for a Ket in
                the same Hilbert space (with the same operator, in this case)
                to perform a|n> = sqrt(n)|n-1> and 'searches backwards' to
                perform <n|a = sqrt(n+1)<n+1|. If it finds another operator,
                it stops."""
                for m in range(n+1, len(term.factors)):
                    if isinstance(term.factors[m], Ket):
                        if term.factors[n].op == term.factors[m].op:
                            term.coef = term.coef*(term.factors[m].eig)**(Add(1)/Add(2))
                            term.factors[m] = Ket(term.factors[m].eig - 1, term.factors[m].op)
                            term.factors[n] = 'out'
                            break
                    elif isinstance(term.factors[m], Operator) or \
                  isinstance(term.factors[m], AdjointOperator) or isinstance(term.factors[m], CreationOperator) or \
                  isinstance(term.factors[m], AnnihilationOperator):
                        break
                if term.factors[n] != 'out':
                    for m in reversed(range(n)):
                        if isinstance(term.factors[m], Bra):
                            if term.factors[n].op == term.factors[m].op:
                                if term.factors[n].maximum == 0:
                                    term.coef = term.coef*(term.factors[m].eig + 1)**(Add(1)/Add(2))
                                    term.factors[m] = Bra(term.factors[m].eig + 1, term.factors[m].op)
                                    term.factors[n] = 'out'
                                else:
                                    if term.factors[m].eig < term.factors[n].maximum:
                                        term.coef = term.coef*(term.factors[m].eig + 1)**(Add(1)/Add(2))
                                        term.factors[m] = Bra(term.factors[m].eig + 1, term.factors[m].op)
                                        term.factors[n] = 'out'
                                    else:
                                        term.coef = Add(0)
                                break
                        elif isinstance(term.factors[m], Operator) or \
                    isinstance(term.factors[m], AdjointOperator) or isinstance(term.factors[m], CreationOperator) or \
                    isinstance(term.factors[m], AnnihilationOperator):
                            break
            else:
                """If a 'factor' inside a term is not a Ket, Bra, operator
                or dual operator, then it is supposed to be a scalar."""
                if term.factors[n] != 'out':
                    term.coef = term.coef*term.factors[n]
                    term.factors[n] = 'out'
        # Remove the 'factors' labeled as 'out':
        while 'out' in term.factors:
            term.factors.remove('out')
    # Arrange each term into the 'canonic form':
    for n in range(len(expression.terms)):
        expression.terms[n] = _canonic_mult(expression.terms[n])
    """Now, when two terms are 'equal', as determined by the function
    is_equal, the first one is replaced by their addition and the second
    one is labeled as 'out', to be removed later. Note that all terms
    of expression belong to the class DiracMult."""
    for n in range(len(expression.terms)):
        for m in range(n+1, len(expression.terms)):
            if _is_equal(expression.terms[n], expression.terms[m]):
                expression.terms[n] = DiracMult(expression.terms[n].coef + \
                                expression.terms[m].coef, *expression.terms[n].factors)
                expression.terms[m] = 'out'
    # Remove terms labeled as 'out':
    while 'out' in expression.terms:
        expression.terms.remove('out')
    # Remove terms with zero coefficient:
    for n in range(len(expression.terms)):
        if expression.terms[n].coef == 0:
            expression.terms[n] = 'out'
    # Remove terms labeled as 'out':
    while 'out' in expression.terms:
        expression.terms.remove('out')
    """Presently, expression consists of a sum of expanded and simplified
    terms, all of them of type DiracMult. Some of them consist of an empty DiracMult
    (empty list) with some non-zero coefficient. These terms are assigned
    the value of their coefficient (they will no longer belong to the DiracMult
    class)."""
    for n in range(len(expression.terms)):
        if len(expression.terms[n].factors) == 0:
            expression.terms[n] = expression.terms[n].coef
    """The expression consists of a sum of DiracMult objects and 'numbers' (i.e.,
    anything which is not a Ket, Bra, Operator, AdjointOperator, DiracMult or DiracSum), and it is
    finally passed to the function _canonic_suma, to order the items in a
    consistent way."""
    expression = _canonic_suma(expression)
    return expression


def _expandir(expression):
    """This function completely expands an expression by repeatedly
    applying the distributive rule. It always returns a sum object
    (DiracSum)."""
    if isinstance(expression, DiracSum):
        nested = True
        while nested:
            n = 0; nmax = len(expression.terms)
            while n <= nmax:
                if n == nmax:
                    nested = False
                    break
                if isinstance(expression.terms[n], DiracMult):
                    k = 0; kmax = len(expression.terms[n].factors)
                    while k < kmax:
                        if isinstance(expression.terms[n].factors[k], DiracSum):
                            L = []
                            for term in expression.terms[n].factors[k].terms:
                                L.append(DiracMult(expression.terms[n].coef, *expression.terms[n].factors[:k] + [term] + expression.terms[n].factors[k+1:]))
                            expression = DiracSum(*expression.terms[:n] + L + expression.terms[n+1:])
                            k = kmax; n = nmax
                        else:
                            k = k + 1
                    n = n + 1
                else:
                    n = n + 1
        return expression
    elif isinstance(expression, DiracMult):
        return _expandir(DiracSum(expression))
    else:
        return DiracSum(expression)


def _canonic_mult(term):
    """This function sorts the elements of a "multiplication" into
    the canonic form, which is: kets followed by bras with their
    corresponding eigen-operators and eigenvalues alphabetically ordered.
    Operators and dual operators remain in place."""
    L = []; kets = []; bras = []
    position_lock = True
# when 'position_lock = False', Kets and Bras can be re-ordered
    for k, u in enumerate(term.factors):
        if (isinstance(u, Operator) or isinstance(u, AdjointOperator) or
isinstance(u, CreationOperator) or isinstance(u, AnnihilationOperator)):
            if not(position_lock):
                kets.sort(); bras.sort(); bras.reverse()
# The lists 'kets' and 'bras' contain the information about the
# accumulated Kets and Bras. When sorted, priority is given to the
# object's eigen-operator or vector space, then to its eigenvalue,
# and finally to its original position in the "multiplication". Of course,
# if a tensor product is to be valid, its kets (or bras) should not live
# in the same space, i.e. they should not have the same eigen-operator.
                for t in kets:
                    # add the ordered kets to the final sequence
                    L.append(term.factors[t[2]])
                for t in bras:
                    # add the ordered bras to the final sequence:
                    L.append(term.factors[t[2]])
                # add the present object (an operator of any kind):
                L.append(u)
                # clear the auxiliary lists and lock the objects' positions:
                kets = []; bras = []
                position_lock = True
            else:
                L.append(u)
        elif isinstance(u, Ket):
            position_lock = False
            kets.append((u.op, u.eig, k))
        elif isinstance(u, Bra):
            position_lock = False
            bras.append((u.op, u.eig, k))
        else:
            print ('Error in canonic_mult. Type not valid.')
    """If the last element of the 'multiplication' is not an operator,
    the final kets and bras must be added:"""
    if len(kets) != 0 or len(bras) != 0:
        kets.sort(); bras.sort(); bras.reverse()
        for t in kets:
            L.append(term.factors[t[2]])
        for t in bras:
            L.append(term.factors[t[2]])

    """The objective of the following lines of code is to prevent
    elements of a tensor product to have the same operator."""
    #for n in range(len(L)-1):
        #if isinstance(L[n], Ket) or isinstance(L[n], Bra):
            #if type(L[n]) == type(L[n+1]) and L[n].op == L[n+1].op:
                #try:
                    #raise Exception, 'Elements of a tensor \
#product cannot be eigenvectors of the same operator.'
                #except Exception:
                    #print 'An exception occurred.'
                    #raise

    term = DiracMult(term.coef, *L)
    return term


def _canonic_suma(expression):
    """This function sorts the terms of a sum into the 'canonic form',
    which is similar to that of the 'multiplication'. The expression
    passed by the function evaluar consists of a sum of non-expandable,
    simplified terms, all of them either 'multiplication' objects (DiracMult),
    or plain, independent terms (numbers; not kets, etc.)."""
    L1 = []; indep = 0
    for k, u in enumerate(expression.terms):
        if isinstance(u, DiracMult):
            if isinstance(u.factors[0], CreationOperator):
                L1.append((1, u.factors[0].op, 0, k))
            if isinstance(u.factors[0], AnnihilationOperator):
                L1.append((2, u.factors[0].op, 0, k))
            if isinstance(u.factors[0], AdjointOperator):
                L1.append((3, u.factors[0].op, 0, k))
            elif isinstance(u.factors[0], Operator):
                L1.append((4, u.factors[0].op, 0, k))
            elif isinstance(u.factors[0], Ket):
                L1.append((5, u.factors[0].op, u.factors[0].eig, k))
            elif isinstance(u.factors[0], Bra):
                L1.append((6, u.factors[0].op, u.factors[0].eig, k))
        else:
            indep = indep + u
    L1.sort(); L2 = []
    if indep != 0:
        L2.append(indep)
    for t in L1:
        L2.append(expression.terms[t[3]])
    if len(L2) == 0:
        expression = 0
    elif len(L2) == 1:
        expression = L2[0]
    else:
        expression = DiracSum(*L2)
    return expression


def _is_equal(u, v):
    """This function returns True if expressions u and v are equal
    and False otherwise. It is supposed that both arguments are
    evaluated terms sorted in the canonical form."""
    if type(u) != type(v):
        return False
    else:
        if u == 'out':
            return False
        elif isinstance(u, Ket) or isinstance(u, Bra):
            if u.op == v.op and u.eig == v.eig:
                return True
            else:
                return False
        elif isinstance(u, Operator) or isinstance(u, AdjointOperator) or isinstance(u, CreationOperator) or isinstance(u, AnnihilationOperator):
            if u.op == v.op:
                return True
            else:
                return False
        elif isinstance(u, DiracMult):
            if len(u.factors) != len(v.factors):
                return False
            else:
                for k in range(len(u.factors)+1):
                    if k == len(u.factors):
                        return True
                        break
                    if _is_equal(u.factors[k], v.factors[k]):
                        continue
                    else:
                        return False
                        break
        else:
            return False