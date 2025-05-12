from sympy import Matrix, I, symbols, Add
from .diracobject import Ket,Bra,Operator,AdjointOperator,_DiracSum,_DiracMult,_expandir

__all__ = ['dag','comm','ptrace','trace','normalize','qmatrix','op2dict','sub_qexpr']

def dag(qobj):
    """
    This function calculates the hermitian conjugate of a quantum object.

    :param qobj: Quantum object. Could be a :obj:`Bra`, :obj:`Ket`, :obj:`Operator` or complex number.
    :type qobj: :obj:`openket.DiracObject`
    :return: The hermitian conjugate of ``qobj``.
    :rtype: :obj:`openket.DiracObject`

    Examples
    ^^^^^^^^^

        .. code-block:: python

            >>> dag(1 + I)
            1-I

        .. code-block:: python

            >>> dag(Ket(0))
            <0|
            >>> dag(Bra(1))
            |1>

    """

    if isinstance(qobj, Bra):
        return Ket(qobj.eig, qobj.op)
    elif isinstance(qobj, Ket):
        return Bra(qobj.eig, qobj.op)
    elif isinstance(qobj, Operator):
        return AdjointOperator(qobj.op)
    elif isinstance(qobj, AdjointOperator):
        return op(qobj.op)
    elif isinstance(qobj, _DiracSum):
        return _DiracSum(*map(dag, qobj.terms))
    elif isinstance(qobj, _DiracMult):
        l = qobj.factors
        L = l[:]
        L.reverse()
        return _DiracMult(qobj.coef.conjugate(), *map(dag, L))
    else:
        return qobj.conjugate()

def comm(A, B):
    """This function calculates the commutator of two operators.

    :param A: Operator one.
    :type A: :obj:`Operator` or matrix
    :param B: Operator two.
    :type B: :obj:`Operator` or matrix
    :return: The conmutator of ``A`` and ``B``.
    :rtype: :obj:`Operator`
    """
    return A*B - B*A

def ptrace(A, htag):
    """This function computes the partial trace of a series of outer products or a bipartite density matrix.

    :param A: Bipartite density matrix or total expression of outer products as linear combination of :obj:`Ket`:obj:`Bra` objects.
    :type A: openket expression
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

def trace(A, basis = 'default'):
    """
    This function computes the total trace of a sum of exterior products or density matrix.
    It finds the total number of Hilbert spaces and their tags and then it uses :obj:`openket.ptrace <openket.core.metrics.ptrace>`
    over all of them succesively, returning a complex number.

    :param A: Bipartite density matrix or total expression of outer products as linear combination of :obj:`Ket`:obj:`Bra` objects.
    :type A: openket expression
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
            A = ptrace(A, i)
    else:
        temp = 0
        for l in basis:
            temp = temp + dag(l)*A*l
        A = temp
    return A

def normalize(state):
    """
    This function produces a vector in a Hilbert space with norm equal to unity.

    :param state: Total expression of vector(s) or a quantum state.
    :type state: :obj:`openket.DiracObject`
    :return: The normalized corresponding state.
    :rtype: :obj:`openket.DiracObject`
    """
    norm = (dag(state)*state)
    norm = float(norm)
    state = state *norm**-0.5
    return state

def qmatrix(A, basis = 'default'):
    """
    This function returns the matrix representation of any operator.

    :param A: Sum of outer products as a linear combination of :obj:`Ket`:obj:`Bra` objects.
    :type A: openket expression
    :param basis: List of :obj:`Ket` objects that conforms the basis in which the matrix will be represented.
                    If no basis is specified, the function founds out how many different kets or bras
                    the expression (sum of the outer products) contains, to determine the basis and the size of the matrix.
                    For the moment, these can only be eigenvectors of a single operator.
    :type basis: array, optional
    :raises Exception: Raised if exterior products are expressed in differents Hilbert spaces.
    :return: The matrix representation of ``A`` in ``basis`` or a default basis.
    :rtype: sympy.Matrix
    """
    A = _DiracSum(A)
    if basis == 'default':
        eigs = []; ops = []
        for term in A.terms:
            if isinstance(term, _DiracMult):
                pass
            else:
                term = _DiracMult(1, term)
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
                row.append(dag(i)*A*j)
            M.append(row)
        M = Matrix(M)
    return M

def op2dict(op, basis):
    """This function creates a dictionary object with keys as :math:`\\langle i|\\rho|j \\rangle` objects in a given basis and their definitions as values.
    The function creates variables representing a complex number, *e.g. a+ib*, where *a* and *b* are defined as real variables.
    If the basis used has :math:`n` elements, :math:`n^2` variables are created.
    It is assumed that all variables are complex and the definition of each :math:`\\langle i|\\rho|j \\rangle` are the variables themselves.

    :param op: The operator involved.
    :type op: :obj:`Operator`
    :param basis: List of :obj:`Ket` objects containing all the elements of the basis.
    :type basis: array
    :return: Dictionary object whose keys correspond to :math:`\\langle i|` A :math:`|j \\rangle` objects in the basis ``basis``,
            and the values to their corresponding defition represented as variable.
    :rtype: dict

    Example
    ^^^^^^^^^

        .. code-block:: python

            >>> b = [Ket(0), Ket(1)]
            >>> A = Operator("A")
            >>> dic = op2dict(op=R, basis=b)
            >>> dic
            {'<0|A|0>': I*Im0 + Re0, '<0|A|1>': I*Im1 + Re1, '<1|A|0>': I*Im2 + Re2, '<1|A|1>': I*Im3 + Re3}

    """
    n = len(basis)
    re = []
    im = []
    t = []
    for i in range( n*n ):
        re.append('Re%d' %i)
        re[i] = symbols('Re%d' %i, real=True, each_char=False)
        im.append('Im%d' %i)
        im[i] = symbols('Im%d' %i, real=True, each_char=False)
    for i in range(len(re)):
        t.append(re[i] + I*im[i])
    D = {}
    k = 0
    for i in range(n):
        for j in range(n):
            D[str(dag(basis[i])*op*basis[j])] = t[k]
            k = k + 1
    return D

def sub_qexpr(qexpr, dic):
    """This function looks within a given expression for keys in the :obj:`openket.op2dict <openket.core.metrics.op2dict>`
    dictionary, and substitutes them with their value.
    It discriminates coefficients and other terms which are not in the dictionary, leaving them the same.

    :param qexpr: Operator expansion in therms of the matrix elements, expressed as sum of :math:`\\langle i|\\rho|j \\rangle`.
    :type qexpr: openket expression
    :param dic: Dictionary object given by the :obj:`openket.op2dict <openket.core.metrics.op2dict>` function.
    :type dic: dict
    :return: The total mathematical expression of ``qexpr`` represented by ``dic`` values.
    :rtype: sympy.expr

    Example
    ^^^^^^^^^

        .. code-block:: python

            >>> b = [Ket(0), Ket(1)]
            >>> A = Operator("A")
            >>> dic = op2dict(op=A, basis=b)
            >>> quantum_expr = Bra(0)*A*Ket(0)+2*I*Bra(0)*A*Ket(1)-3.4*Bra(1)*A*Ket(1)
            >>> symbolic_expr = sub_qexpr(qexpr=quantum_expr, dic=dic)
            >>> symbolic_expr
            I*Im0 + Re0 + 2*I*(I*Im1 + Re1) - 3.4*I*Im3 - 3.4*Re3

    """
    if qexpr == 0.0:
        l = [0]
    else:
        qexpr = _DiracSum(qexpr)
        L = qexpr.terms
        l = []
        for i in L:
            y = _DiracMult(1,*i.factors)
            if str(y) in dic:
                y = str(y)
                y = dic[y]
            l.append(y * (i.coef))
    return Add(*l)
