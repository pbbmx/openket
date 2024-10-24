from sympy import Matrix, I, symbols
from openket.core.diracobject import _expandir
from openket.core.diracobject import *
from openket.core.gates import *

def Adj(A):
    """
    This function calculates the hermitian conjugate of a quantum object.

    :param A: Openket expression, this could be a :obj:`Bra`, :obj:`Ket`, :obj:`Operator` or complex number.
    :type A: :obj:`DiracObject`
    :return: The hermitian conjugate of ``A``.
    :rtype: :obj:`DiracObject`

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

    Of course, it's all possible to calculate the norm of any state.

        .. code-block:: python

            >>> w = I*Ket(0) - 4*Ket(1)
            >>> w
            I|0> - 4|1>
            >>> norm = Adj(w)*w
            >>> norm
            17

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

    :param A: The operator or matrix.
    :type A: :obj:`Operator`
    :param B: The operator or matrix.
    :type B: :obj:`Operator`
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

    Example
    ^^^^^^^^^

        .. code-block:: python

            >>> psi = Ket(0,"A")*Ket(1,"B") - Ket(1,"A")*Ket(0,"B")
            >>> psi
            |0_A>|1_B> - |1_A>|0_B>
            >>> R = psi*Adj(psi)
            >>> R
            |0_A>|1_B><1_B|<0_A| - |0_A>|1_B><0_B|<1_A| - |1_A>|0_B><1_B|<0_A| + |1_A>|0_B><0_B|<1_A|
            >>>
            >>>
            >>> RB = TraceOut(R,"A")
            >>> RB
            |0_B><0_B| + |1_B><1_B|
            >>> RA = TraceOut(R,"B")
            >>> RA
            |0_A><0_A| + |1_A><1_A|

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
    It finds the total number of Hilbert spaces and their tags and then it uses :obj:`TraceOut` over
    all of them succesively, returning a complex number.

    :param A: Bipartite density matrix or total expression of outer products.
    :type A: Linear combination of :obj:`Ket`:obj:`Bra` objects
    :param basis: 
    :type basis:
    :return: The total trace of ``A`` acting in all Hilbert spaces.
    :rtype: number

    Example
    ^^^^^^^^^
    Continuing with the previus example:

        .. code-block:: python

            >>> Trace(R)
            2

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
    :type state: :obj:`DiracObject`
    :return: The normalized corresponding state.
    :rtype: :obj:`DiracObject`
    """
    norm = (Adj(state)*state)
    norm = float(norm)
    state = state / (norm**0.5)
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
    :type basis: list, optional
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