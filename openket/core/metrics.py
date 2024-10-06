from sympy import Matrix, I, symbols
import openket.core.diracobject as do

def Adj(A):
    """
    This function calculates the hermitian conjugate of a quantum object.

    :param A: Openket expression, this could be a Bra, Ket, Operator or complex number.
    :type A: DiracObject
    :return: The hermitian conjugate of `A`.
    :rtype: DiracObject
    """

    if isinstance(A, do.Bra):
        return do.Ket(A.eig, A.op)
    elif isinstance(A, do.Ket):
        return do.Bra(A.eig, A.op)
    elif isinstance(A, do.Operator):
        return do.AdjointOperator(A.op)
    elif isinstance(A, do.AdjointOperator):
        return do.op(A.op)
    elif isinstance(A, do.DiracSum):
        return do.DiracSum(*map(Adj, A.terms))
    elif isinstance(A, do.DiracMult):
        l = A.factors
        L = l[:]
        L.reverse()
        return do.DiracMult(A.coef.conjugate(), *map(Adj, L))
    else:
        return A.conjugate()

def Commutator(A, B):
    """This function calculates the conmutator of two operators.

    :param A: Operator or matrix.
    :type A: Operator
    :param B: Operator or matrix.
    :type B: Operator
    :return: The conmutator of `A` and `B`.
    :rtype: Operator
    """
    return A*B - B*A

def TraceOut(A, htag):
    """This function computes the partial trace of a series of outer products or a bipartite density matrix.

    :param A: Bipartite density matrix or total expression of outer products.
    :type A: Operator
    :param htag: Tag of the Hilbert space you will be tracing out.
    :type htag: string
    :return: The partial trace of `A` acting only in `htag` space.
    :rtype: Operator
    """
    A = do._expandir(A)
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
                term = term + do.Bra(j,i) * A * do.Ket(j,i)
            A = term
            term = 0
    return A

def Trace(A, basis = 'default'):
    """
    This function computes the total trace of a sum of exterior products or density matrix.
    It finds the total number of Hilbert spaces and their tags and then it uses TraceOut over
    all of them succesively, returning a complex number.

    :param A: Bipartite density matrix or total expression of outer products.
    :type A: Operator
    :param basis: 
    :type basis:
    :return: The total trace of `A` acting in all Hilbert spaces.
    :rtype: number
    """
    if basis == 'default':
        A = do._expandir(A)
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
    :type state: DiracObject
    :return: The normalized corresponding state.
    :rtype: DiracObject
    """
    norm = (Adj(state)*state)
    norm = float(norm)
    state = state / (norm**0.5)
    return state

def Qmatrix(A, basis = 'default'):
    """
    This function returns the matrix representation of any operator.

    :param A: Sum of exterior products.
    :type A: DiracObject
    :param basis: List of Ket objects that conforms the basis in which the matrix will be represented.
                    If no basis is specified, the function founds out how many different kets or bras
                    the expression (sum of the exterior products) contains, to determine the basis and the size of the matrix.
                    For the moment, these can only be eigenvectors of a single operator.
    :type basis: list, optional
    :raises Exception: Raised if exterior products are expressed in differents Hilbert spaces.
    :return: The matrix representation of A in `basis` or a default basis.
    :rtype: Sympy Matrix
    """
    A = do.DiracSum(A)
    if basis == 'default':
        eigs = []; ops = []
        for term in A.terms:
            if isinstance(term, do.DiracMult):
                pass
            else:
                term = do.DiracMult(1, term)
            for fact in term.factors:
                if isinstance(fact, do.Ket) or isinstance(fact, do.Bra):
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
                row.append(do.Bra(i,oper)*A*do.Ket(j,oper))
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
    """This function creates a dictionary object with 'words' as strings of
    matrix elements of the operator R, in the base x. x is given as a list
    containing all the elements of the basis, and R is necessarily defined
    as an Operator. The function creates variables representing a complex
    number, e.g. w + I*y, where w and y are defined as real variables. If the
    base used has n elements, n*n variables are created, asigning each
    'word' a variable. Is is assumed that all variables al complex.
    The 'definitions' of each 'word' are the variables themselves."""
    """
    This function creates a dictionary object with 'words' as strings of matrix elements of the operator `A`, in the base `basis`.

    :param A: _description_
    :type A: Operator
    :param basis: List of Ket objects containing all the elements of the basis.
    :type basis: list
    :return: _description_
    :rtype: _type_
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