from sympy import Matrix, I, symbols, Add
from diracobject import AdjointOperator, DiracSum, DiracMult, _evaluar, _expandir, _canonic_mult, _canonic_suma, _is_equal
from diracobject import *
from gates import *

def Adj(A):
    """
    This function calculates the hermitian conjugate of a quantum object.

    :param A: Could be a :obj:`Bra`, :obj:`Ket`, :obj:`Operator` or complex number.
    :type A: :obj:`openket.DiracObject`
    :return: The hermitian conjugate of ``A``.
    :rtype: :obj:`openket.DiracObject`

    Examples
    ^^^^^^^^^

        .. code-block:: python

            >>> Adj(1 + I)
            1-I

        .. code-block:: python

            >>> Adj(Ket(0))
            <0|
            >>> Adj(Bra(1))
            |1>

    """

    if isinstance(A, Bra):
        return Ket(A.eig, A.op)
    elif isinstance(A, Ket):
        return Bra(A.eig, A.op)
    elif isinstance(A, Operator):
        return AdjointOperator(A.op)
    elif isinstance(A, AdjointOperator):
        return op(A.op)
    elif isinstance(A, DiracSum):
        return DiracSum(*map(Adj, A.terms))
    elif isinstance(A, DiracMult):
        l = A.factors
        L = l[:]
        L.reverse()
        return DiracMult(A.coef.conjugate(), *map(Adj, L))
    else:
        return A.conjugate()

def Commutator(A, B):
    """This function calculates the conmutator of two operators.

    :param A: Operator one.
    :type A: :obj:`Operator` or matrix
    :param B: Operator two.
    :type B: :obj:`Operator` or matrix
    :return: The conmutator of ``A`` and ``B``.
    :rtype: :obj:`Operator`
    """
    return A*B - B*A

def TraceOut(A, htag):
    """This function computes the partial trace of a series of outer products or a bipartite density matrix.

    :param A: Bipartite density matrix or total expression of outer products.
    :type A: Linear combination of :obj:`Ket`:obj:`Bra` objects
    :param htag: Tag of the Hilbert space you will be tracing out.
    :type htag: string
    :return: The partial trace of ``A`` acting only in ``htag`` space.
    :rtype: Linear combination of :obj:`Ket`:obj:`Bra` objects

    """
    A = _expandir(A)
    Lterms = A.terms
    htags = []
    elements = []
    l = []
    for i in Lterms:
        Lfactors = i.factors
        for j in Lfactors:
            htags.append(j.op)
            elements.append(j)
    htags = list(set(htags))
    d = {}
    for i in htags:
        d[i] = []
    for i in d:
        for j in elements:
            if j.op == i:
                d[i].append(j.eig)
    for i in d:
        d[i] = list(set(d[i]))
    term = 0
    for i in d.keys():
        if i == htag:
            for j in d[i]:
                term = term + Bra(j,i) * A * Ket(j,i)
            A = term
            term = 0
    return A

def Trace(A, basis = 'default'):
    """
    This function computes the total trace of a sum of exterior products or density matrix.
    It finds the total number of Hilbert spaces and their tags and then it uses :obj:`openket.TraceOut <openket.core.metrics.TraceOut>`
    over all of them succesively, returning a complex number.

    :param A: Bipartite density matrix or total expression of outer products.
    :type A: Linear combination of :obj:`Ket`:obj:`Bra` objects
    :param basis: 
    :type basis:
    :return: The total trace of ``A`` acting in all Hilbert spaces.
    :rtype: float

    """
    if basis == 'default':
        A = _expandir(A)
        Lterms = A.terms
        htags = []
        for i in Lterms:
            Lfactors = i.factors
            for j in Lfactors:
                htags.append(j.op)
        htags = list(set(htags))
        for i in htags:
            A = TraceOut(A, i)
    else:
        temp = 0
        for l in basis:
            temp = temp + Adj(l)*A*l
        A = temp
    return A

def Normalize(state):
    """
    This function produces a vector in a Hilbert space with norm equal to unity.

    :param state: Total expression of vector(s) or a quantum state.
    :type state: :obj:`openket.DiracObject`
    :return: The normalized corresponding state.
    :rtype: :obj:`openket.DiracObject`
    """
    norm = (Adj(state)*state)
    norm = float(norm)
    state = state *norm**-0.5
    return state

def Qmatrix(A, basis = 'default'):
    """
    This function returns the matrix representation of any operator.

    :param A: Sum of outer products.
    :type A: Linear combination of :obj:`Ket`:obj:`Bra` objects
    :param basis: List of :obj:`Ket` objects that conforms the basis in which the matrix will be represented.
                    If no basis is specified, the function founds out how many different kets or bras
                    the expression (sum of the outer products) contains, to determine the basis and the size of the matrix.
                    For the moment, these can only be eigenvectors of a single operator.
    :type basis: array, optional
    :raises Exception: Raised if exterior products are expressed in differents Hilbert spaces.
    :return: The matrix representation of ``A`` in ``basis`` or a default basis.
    :rtype: Sympy Matrix
    """
    A = DiracSum(A)
    if basis == 'default':
        eigs = []; ops = []
        for term in A.terms:
            if isinstance(term, DiracMult):
                pass
            else:
                term = DiracMult(1, term)
            for fact in term.factors:
                if isinstance(fact, Ket) or isinstance(fact, Bra):
                    eigs.append(fact.eig)
                    ops.append(fact.op)
        eigs = list(set(eigs)); eigs.sort()
        ops = list(set(ops)); ops.sort()
        if len(ops) > 1:
            try:
                raise Exception('Matrix form can only be obtained for operators acting on a single Hilbert space.')
            except Exception:
                print('An exception occurred.')
                raise
        else:
            oper = ops[0]
        M = []
        for i in eigs:
            row = []
            for j in eigs:
                row.append(Bra(i,oper)*A*Ket(j,oper))
            M.append(row)
        M = Matrix(M)
    else:
        M = []
        for i in basis:
            row = []
            for j in basis:
                row.append(Adj(i)*A*j)
            M.append(row)
        M = Matrix(M)
    return M

def Dictionary(A, basis):
    """This function creates a dictionary object with keys as :math:`\\langle i|\\rho|j \\rangle` objects in a given basis and their definitions as values.
    The function creates variables representing a complex number, *e.g. a+ib*, where *a* and *b* are defined as real variables.
    If the basis used has :math:`n` elements, :math:`n^2` variables are created.
    It is assumed that all variables are complex and the definition of each :math:`\\langle i|\\rho|j \\rangle` are the variables themselves.

    :param A: The operator involved.
    :type A: :obj:`Operator`
    :param basis: List of :obj:`Ket` objects containing all the elements of the basis.
    :type basis: array
    :return: Dictionary object whose keys correspond to :math:`\\langle i|` ``A`` :math:`|j \\rangle` objects in the basis ``basis``,
            and the values to their corresponding defition represented as variable.
    :rtype: dict

    Example
    ^^^^^^^^^

        .. code-block:: python

            >>> b = [Ket(0), Ket(1)]
            >>> R = Operator("R")
            >>> Dictionary(R, b)
            {'<0|R|0>': I*Im0 + Re0, '<0|R|1>': I*Im1 + Re1, '<1|R|0>': I*Im2 + Re2, '<1|R|1>': I*Im3 + Re3}

    """
    n = len(basis)
    u = []
    w = []
    t = []
    for i in range( n*n ):
        u.append('Re%d' %i)
        u[i] = symbols('Re%d' %i, real=True, each_char=False)
        w.append('Im%d' %i)
        w[i] = symbols('Im%d' %i, real=True, each_char=False)
    for i in range(len(u)):
        t.append(u[i] + I*w[i])
    D = {}
    k = 0
    for i in range(n):
        for j in range(n):
            D[str(Adj(basis[i])*A*basis[j])] = t[k]
            k = k + 1
    return D

def Qch(expr, dic):
    """This function looks within a given expression for keys in the :obj:`openket.Dictionary <openket.core.metrics.Dictionary>`
    dictionary, and substitutes them with their value.
    It discriminates coefficients and other terms which are not in the dictionary, leaving them the same.

    :param expr: Openket expression of a sum of :math:`\\langle i|\\rho|j \\rangle` objects.
    :type expr: Lineal combination of :obj:`Bra`:obj:`Operator`:obj:`Ket` objects.
    :param dic: Dictionary object given by the :obj:`openket.Dictionary <openket.core.metrics.Dictionary>` function.
    :type dic: dict
    :return: The total mathematical expression of ``expr`` represented by ``D`` values.
    :rtype: Sympy expression

    Example
    ^^^^^^^^^

        .. code-block:: python

            >>> b = [Ket(0), Ket(1)]
            >>> R = Operator("R")
            #{'<0|R|0>': I*Im0 + Re0, '<0|R|1>': I*Im1 + Re1, '<1|R|0>': I*Im2 + Re2, '<1|R|1>': I*Im3 + Re3}
            >>> D = Dictionary(R, b)
            #<0|R|0> + 2*I<0|R|1> +  - 3.4<1|R|1>
            >>> expression =  Bra(0)*R*Ket(0)+2*I*Bra(0)*R*Ket(1)-3.4*Bra(1)*R*Ket(1)
            >>> Qch(expr, D)
            I*Im0 + Re0 + 2*I*(I*Im1 + Re1) - 3.4*I*Im3 - 3.4*Re3

    """
    if expr == 0.0:
        l = [0]
    else:
        expr = DiracSum(expr)
        L = expr.terms
        l = []
        for i in L:
            y = DiracMult(1,*i.factors)
            if str(y) in dic:
                y = str(y)
                y = dic[y]
            l.append(y * (i.coef))
    return Add(*l)
