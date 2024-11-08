from openket.core.diracobject import *
from openket.core.metrics import *
from sympy import Symbol, symbols, expand
from sympy.core.cache import clear_cache
from pylab import conjugate
import re as regex


def Dictionary(A, basis):
    """This function creates a dictionary object with keys as :math:`\\langle i|\\rho|j \\rangle` objects in a given basis and their definitions as values.
    The function creates variables representing a complex number, *e.g. a+ib*, where *a* and *b* are defined as real variables.
    If the basis used has :math:`n` elements, :math:`n^2` variables are created.
    It is assumed that all variables are complex and the definition of each :math:`\\langle i|\\rho|j \\rangle` are the variables themselves.

    :param A: The operator involved.
    :type A: :obj:`Operator`
    :param basis: List of :obj:`Ket` objects containing all the elements of the basis.
    :type basis: list of :obj:`Ket`
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
    """This function looks within a given expression for keys in the :obj:`Dictionary` dictionary, and substitutes them with their value.
    It discriminates coefficients and other terms which are not in the dictionary, leaving them the same.

    :param expr: Openket expression of a sum of :math:`\\langle i|\\rho|j \\rangle` objects.
    :type expr: Lineal combination of :obj:`Bra`:obj:`Operator`:obj:`Ket` objects.
    :param dic: Dictionary object given by the :obj:`Dictionary` function.
    :type dic: dict
    :return: The total mathematical expression of ``expr`` represented by ``D`` values.
    :rtype: Mathematical expression with real variables

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

def Qeq(R, Rdot, basis, filename = "func", dictname = "dic", GSL = False):

    """
    The main purpose for this function is to study time evolution of an :obj:`Operator`, with an expression
    for its time derivative, both given in a specified basis.
    This function creates a file called ``filename``, which defines a function
    in the sintaxis needed to feed Scipy's ODE solver `odeint` or GSL ODE solver.

    :param R: Operator involved. The function turns it in terms of the basis ``basis``.
    :type R: :obj:`Operator`
    :param Rdot: Operator that expresses change over the time, usually the commutator [H,R] where H is the Hamiltonian of the system.
    :type Rdot: :obj:`Operator`
    :param basis: List of :obj:`Ket` objects containing all the elements of the basis.
    :type basis: list of :obj:`Ket`
    :param filename: The name of the out file with the function.
    :type filename: string
    :param dictname: The name of the dictionay file.
    :type dictname: string
    :param GSL: Defines the syntaxis of the ODE solver; if :code:`GLS = True` define a function in the sintaxis for GSL ODE solver,
                if :code:`False` define in the sintaxis needed to feed Scipy's ODE solver `odeint`.
    :type GSL: bool
    :return: Return a file with a function written in the sintaxis defined. Besides, if :code:`GLS = True` also returns a
            file with the dictionary of the base and the variables.
            It first creates a dictionary with the operator ``R`` elements expressed in the basis ``basis``
            using the function :obj:`Dictionary`, and then it finds the different expressions of time evolution 
            for every matrix element substituting them for variables using :obj:`Qch`.
            Finally it writes the function :math:`f` into a file, where the ODE has the form :math:`\\frac{dy}{dt} = f(y,t)`;
            :math:`y` should be understood as a vector, which is defined entry-like as a list,and :math:`t` is the time variable.
            :math:`f` returns also a list, containig the expressions of the derivatives of each entry
            of :math:`y` as entries of another list. 
            Scipy's ODE solver `odeint` understands this as a collection of coupled ODE to be solved, specifying
            an initial condition and the time step.
    :rtype: file

    Example
    ^^^^^^^^^

    We will considerate the Two Level Semiclassical Atoms problem, when an atom interacts with an electromagnetic field.
    We can describe the atom using a 2 element base, let this basis be {:math:`|0 \\rangle, |1 \\rangle`}.
    Furthermore let the energy of the levels be :math:`\\hbar \\omega_0` and :math:`\\hbar \\omega_1` respectively.
    Now, the total Hamiltonian may be divided into two parts: the Hamiltonian of the unperturbed atom :math:`\\hat{H_0}`,
    and the interation part :math:`\\hat{V}`. These can be found to have the following form
    
        .. math::

            \\hat{H_0} = \\hbar \\omega_0 |0 \\rangle \\langle 0| + \\hbar \\omega_1 |1 \\rangle \\langle 1| \\\\
            \\hat{V} = g E(t) |0 \\rangle \\langle 1| + g* E(t) |1 \\rangle \\langle 0|
    
    Where we have put :math:`g = e \\langle 0| \\vec{r} |1 \\rangle \\cdot \\vec{\\epsilon}` and 
    :math:`\\vec{E}(t) = |\\vec{E}(t)|\\vec{\\epsilon} = E(t) \\vec{\\epsilon}` with :math:`|\\vec{\\epsilon}| = 1`.
    We assume the field has the form

        .. math::

            E(t) = A cos(\\nu t)

    It is convenient to change to the interaction picture, and following the rotating wave approximation (RWA)
    we end up with the Hamiltonian

        .. math::

            \\hat{H} = \\hbar (\\Omega |0 \\rangle \\langle 1| + \\Omega * |1 \\rangle \\langle 0|)

    Where :math:`\\Omega = \\frac{gA}{2 \\hbar}` and we assume :math:`\\Omega = 1` for simplicity.
    Following is the script which numerically solves the time evolution of the system described.

        .. code-block:: python
        
            >>> import matplotlib.pyplot as plt
            >>> from numpy import linspace
            >>> from scipy.integrate import odeint

            >>> basis = [Ket(0), Ket(1)]
            >>> rho = Operator("R")
            >>> rho_dot = -I*Commutator(H, rho)
            >>> H = Ket(0)*Bra(1) + Ket(1)*Bra(0)
            >>> Qeq(rho, rho_dot, basis, "func", "dic")

            # The function created a file called `func` with the function we have put in `odeint`
            >>> with open("func") as file:
                    exec(file.read())
            
            >>> initial_conditions = [1,0, 0,0, 0,0, 0,0] #[1+i0, 0+i0, 0+i0, 0+i0]
            >>> t = linspace(0,10,1000)
            >>> solution = odeint(f, initial_conditions, t)
            >>> plt.plot(t, solution[:,3])
            >>> plt.show()

    Note that the initial condition is given as a list, where the position in the list indicates the bracket term that one wishes
    to give as the initial condition (the dic tells you which is which). As we are working with complex numbers, we represent one number
    with two positions in the list.
    Finally, the command :code:`solution[:,3]` partitions the total solution and only keeps the third entry 
    (:code:`y[3]` corresponding to :math:`\\langle 1| \\rho |1 \\rangle`) in all the lists within `solution`.
    The graph we obtain is the following.

    .. figure:: qeq-example.png
        :scale: 85 %

        Time evolution for :math:`\\langle 1| \\rho |1 \\rangle` with the initial condition :math:`\\rho_{t=0} = |0 \\rangle \\langle 0|`
    
    """
    
    #GSL ODE solver (C language)
    if GSL:  
        #se crea el diccionario para la base b y el operador R
        D = Dictionary(R, basis)
        f = open(filename, 'w')
        f.write('#include <stdio.h>\n') 
        f.write('#include <complex.h>\n') 
        f.write('#include <gsl/gsl_errno.h>\n')
        f.write('#include <gsl/gsl_matrix.h>\n')
        f.write('#include <gsl/gsl_odeiv2.h>\n\n')
        f.write('int func (double t, double y[], double f[], void *params) {\n')

        # crea la lista de valores reales e imaginarios W = [Re1,Im1,...,...,Ren*n,Imn*n]
        W = []
        n = len(basis)
        for i in range( n*n ):
            W.append(symbols('Re%d' %k, real=True, each_char=False))
            W.append(symbols('Im%d' %k, real=True, each_char=False))

        # crea una lista con los strings L = [y[0],y[1],...,y[2*n*n]]
        L = ['y[' + '%d' %i + ']' for i in range(2*n*n)]

        # crea el diccionario lista = {'Re1': y[0], ... ,'Imn*n': y[2n*n]}
        lista={}
        for mi in range(2*n*n): lista[W[mi]]=Symbol(L[mi])

        # escribe el sistema de ecuaciones para encontrar los valores a,b de los coeficientes a+ib
        for i,ele in enumerate(basis):
            bra = Adj(ele)
            print (i,"/",len(basis))
            clear_cache()
            for j,ket in enumerate(basis):
                k = 2*(i*n + j)
                z = bra*Rdot*ket
                zz = Qch(z, D)
                re = expand((zz+conjugate(zz))/2.0,complex=True)
                im = expand((zz-conjugate(zz))/(2.0*I),complex=True)
                for var_subs in re.atoms(Symbol):
                    if var_subs in lista:
                        re = re.subs(var_subs,lista[var_subs])
                f.write('  ' + 'f[%d] = ' %k + str(re) + ';\n')
                for var_subs in im.atoms(Symbol):
                    if "var_subs" in lista:
                        im = im.subs(var_subs,lista[var_subs])
                k = k + 1
                f.write('  ' + 'f[%d] = ' %k + str(im) + ';\n')
        f.write('    ' + 'return GSL_SUCCESS;\n')
        f.write('}\n\n')
        dim = len(L)
        f.write('int dim = ' + str(dim) + ';\n')
        f.close()

        # crea el diccionario que relaciona '<i|R|j>' con las y[i]
        f = open(dictname, 'w')
        f.write(dictname + '=' + '{')
        for bra in basis:
            bra = Adj(bra)
            for ket in basis:
                temp = bra*R*ket
                temp = Qch(temp, D)
                for var_subs in temp.atoms(Symbol):
                    etiqueta=regex.search("\[(.*)\]",str(lista[var_subs]))
                    temp = temp.subs(var_subs,Symbol("var('y"+etiqueta.groups()[0]+"',real=True)"))
                f.write("'" + str(bra*R*ket) + "'" + ':' + \
                str(temp) + ',\n')
        f.write( '}')
        f.close()

    #Scipy's ODE solver (python language)
    else:
        D = Dictionary(R, basis)
        f = open(filename, 'w')
        f.write('from sympy import var, I\n')
        f.write('def f(y, t):\n')
        
        # crea la lista de valores reales e imaginarios W = [Re1,Im1,...,...,Ren*n,Imn*n]
        W = []
        n = len(basis)
        for k in range(n*n):
            W.append(symbols('Re%d' %k, real=True, each_char=False))
            W.append(symbols('Im%d' %k, real=True, each_char=False))

        # crea una lista con los strings L = [y[0],y[1],...,y[2*n*n]]
        L = ['y[' + '%d' %i + ']' for i in range(2*n*n)]

        # crea el diccionario lista = {'Re1': y[0], ... ,'Imn*n': y[2n*n]}
        lista={}
        for mi in range(2*n*n): lista[W[mi]]=Symbol(L[mi])
        
        # escribe el sistema de ecuaciones para encontrar los valores a,b de los coeficientes a+ib
        f.write('    ' + 'return [')    
        for i,ele in enumerate(basis):
            bra = Adj(ele)
            print (i,"/",len(basis))
            clear_cache()
            for j,ket in enumerate(basis):
                z = bra*Rdot*ket
                zz = Qch(z, D)
                re = expand((zz+conjugate(zz))/2.0,complex=True)
                im = expand((zz-conjugate(zz))/(2.0*I),complex=True)
                for var_subs in re.atoms(Symbol):
                    if var_subs in lista:
                        re = re.subs(var_subs,lista[var_subs])
                f.write(str(re)+',\n')
                for var_subs in im.atoms(Symbol):
                    if var_subs in lista:
                        im = im.subs(var_subs,lista[var_subs])
                f.write(str(im))
                if basis.index(ele)*basis.index(ket)==(len(basis)-1)**2:
                    f.write(']\n')
                else:
                    f.write(',\n')
        f.close()

        # crea el diccionario que relaciona '<i|R|j>' con las y[i]
        f = open(dictname, 'w')
        for bra in basis:
            bra = Adj(bra)
            for ket in basis:
                temp = bra*R*ket
                temp = Qch(temp, D)
                for var_subs in temp.atoms(Symbol):
                    etiqueta = regex.search("\[(.*)\]",str(lista[var_subs]))
                    temp = temp.subs(var_subs,Symbol("var('y"+etiqueta.groups()[0]+"',real=True)"))
                f.write("'" + str(bra*R*ket) + "'" + ':' + \
                str(temp) + ',')
        f.write( '}')
        f.close()
