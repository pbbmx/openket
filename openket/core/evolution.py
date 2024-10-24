from openket.core.diracobject import *
from openket.core.metrics import *
from sympy import Symbol, symbols, expand
from sympy.core.cache import clear_cache
from pylab import conjugate
import re


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

    :param R: Operator involved. It should be in terms of the basis ``basis``.
    :type R: :obj:`Operator`
    :param Rdot: Usually the commutator of R with H, where H is the Hamiltonian of the system.
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
            It first creates a dictionary with the function :obj:`Dictionary`, using the list ``basis``,
            and the operator ``R``, and then it finds the different expressions of time evolution 
            for every matrix element substituting them for variables using :obj:`Qch`.
            Finally it writes the function :math:`f` into a file, where the ODE has the form :math:`\\frac{dy}{dt} = f(y,t)`;
            :math:`y` should be understood as a vector, which is defined entry-like as a list, for example :math:`y[0]=x, y[1]=w`, etc.,
            and :math:`t` is the time variable. :math:`f` returns also a list, containig the expressions of the derivatives of each entry
            of :math:`y` as entries of another list. 
            Scipy's ODE solver `odeint` understands this as a collection of coupled ODE to be solved, specifying
            an initial condition and the time step.
    :rtype: file
    """
    
    #GSL ODE solver
    if GSL:  
        #se crea el diccionario para la base b y el operador R
        D = Dictionary(R, basis)
        #se empieza a escribir la funcion para lenguaje C
        f = open(filename, 'w')
        f.write('#include <stdio.h>\n') 
        f.write('#include <complex.h>\n') 
        f.write('#include <gsl/gsl_errno.h>\n')
        f.write('#include <gsl/gsl_matrix.h>\n')
        f.write('#include <gsl/gsl_odeiv2.h>\n\n')
        # f.write('#include "hdf5.h"\n')
        # f.write('#define FILE "dset.h5"\n\n')

        f.write('int func (double t, double y[], double f[], void *params) {\n')

        #se crea la lista de valores reales e imaginarios W = [Re1,Im1,...,...,Ren*n,Imn*n]

        W = []
        a = []
        d = []
        n = len(basis)

        for i in range( n*n ):
            a.append('Re%d' %i)
            a[i] = symbols('Re%d' %i, real=True, each_char=False)
            d.append('Im%d' %i)
            d[i] = symbols('Im%d' %i, real=True, each_char=False)
        k = 0
        for i in range(n):
            for j in range(n):
                W.append(a[k])
                W.append(d[k])
                k = k + 1

        #termina de crearse W = [Re1,Im1,...,...,Ren*n,Imn*n]

        L = ['y[' + '%d' %i + ']' for i in range(2*n*n)]

        # Crea L = [y[1],...,y[2*n*n]]

        # Crea el diccionario para hacer la substitucion de forma rapida
        lista={}
        for mi in range(2*n*n):
            lista[W[mi]]=Symbol(L[mi])

        # Termina de crear un diccionario, lista = {'Re1': y[0], ... ,'Imn*n': y[2n*n]}  

        for numero1,i in enumerate(basis): #numero es el numero del elemento, i es el elemento i-esimo en b
            print (numero1,"/",len(basis))
            clear_cache()   #Limpia el cache para que no ocupe toda la memoria
            for numero2,j in enumerate(basis):
                k =2*(numero1*(len(basis)) + numero2)
                z = Adj(i)*Rdot*j
                zz = Qch(z, D)
                zz = expand((zz+conjugate(zz))/2.0,complex=True) #Parte real
                for var_subs in zz.atoms(Symbol):
                    if var_subs in lista:
                        zz = zz.subs(var_subs,lista[var_subs])
                f.write('  ' + 'f[%d] = ' %k + str(zz) + ';\n')         
                zz = Qch(z, D)
                zz = expand((zz-conjugate(zz))/(2.0*I),complex=True) #Parte imaginaria
                for var_subs in zz.atoms(Symbol):
                    if lista.has_key(var_subs):
                        zz = zz.subs(var_subs,lista[var_subs])
                k = k + 1
                f.write('  ' + 'f[%d] = ' %k + str(zz) + ';\n')     
                

        f.write('    ' + 'return GSL_SUCCESS;\n')  
        f.write('}\n\n')
        dim = len(L)
        f.write('int dim = ' + str(dim) + ';\n') 
        f.close()  

        # Se termina de escribir la funcion  


        # el diccionario que relaciona '<i|R|j>' con las y[i]
        # realmente con las "var('y[i]')"

        f = open(dictname, 'w')
        f.write(dictname + '=' + '{')
        for count1 in basis:
            for count2 in basis:
                temp = Adj(count1)*R*count2
                temp = Qch(temp, D)
                for var_subs in temp.atoms(Symbol):
                    etiqueta=re.search("\[(.*)\]",str(lista[var_subs]))
                    temp = temp.subs(var_subs,Symbol("var('y"+etiqueta.groups()[0]+"',real=True)"))
                f.write("'" + str(Adj(count1)*R*count2) + "'" + ':' + \
                str(temp) + ',\n')
        f.write( '}')

        f.close()
    #Scipy's ODE solver
    else:
    
        D = Dictionary(basis, R)
        #print D
        f = open(filename, 'w')
        f.write('from sympy import var, I\n')
        f.write('def f(y, t):\n')
        
        # crea la lista con las variables Re, Im
        M, u, w = [], [], []
        n = len(basis)
        for i in range( n*n ):
            u.append('Re%d' %i)
            u[i] = symbols('Re%d' %i, real=True, each_char=False)
            w.append('Im%d' %i)
            w[i] = symbols('Im%d' %i, real=True, each_char=False)
        k = 0
        for i in range(n):
            for j in range(n):
                M.append(u[k])
                M.append(w[k])
                k = k + 1
        # termina de crearla

        # crea la lista con los strings y[0], y[1]
        L = ['y[' + '%d' %i + ']' for i in range(2*n*n)]
        # termina de crearla

        # Crea el diccionario para hacer la substitucion de forma rapida
        lista={}
        for mi in range(2*n*n):
            lista[M[mi]]=Symbol(L[mi])

        f.write('    ' + 'return [')    
        for numero,i in enumerate(basis):
            print (numero,"/",len(basis))
            clear_cache()   #Limpia el cache para que no ocupe toda la memoria
            for j in basis:
                z = Adj(i)*Rdot*j
                zz = Qch(z, D)
                zz = expand((zz+conjugate(zz))/2.0,complex=True)
                for var_subs in zz.atoms(Symbol):
                    if var_subs in lista:
                        zz = zz.subs(var_subs,lista[var_subs])
                f.write( str(zz)+',\n')         
                zz = Qch(z, D)
                zz = expand((zz-conjugate(zz))/(2.0*I),complex=True)
                for var_subs in zz.atoms(Symbol):
                    if var_subs in lista:
                        zz = zz.subs(var_subs,lista[var_subs])
                f.write(str(zz))     
                if basis.index(i)*basis.index(j)==(len(basis)-1)**2:
                    f.write(']\n')
                else:
                    f.write(',\n')
        # el diccionario que relaciona '<i|R|j>' con las y[i]
        # realmente con las "var('y[i]')"

    ####    EMPIEZA VIEJO
        #for conta1 in D.keys():
        #    for conta2 in range(n**2):
        #        D[conta1] = D[conta1].subs(M[conta2], Symbol(L[conta2])) #cambie a LL
        #for conta3 in D.keys():
        #    D[conta3] = str(D[conta3])
        # termina de crearlo
        #print D
    ####    TERMINA VIEJO    

    #    LL = ["var('y" + '%d' %i + "',real=True)" for i in range(n**2)]
        f.write(dictname + '=' + '{')
        for count1 in basis:
            for count2 in basis:
                temp = Adj(count1)*R*count2
                temp = Qch(temp, D)
                for var_subs in temp.atoms(Symbol):
                    etiqueta=re.search("\[(.*)\]",str(lista[var_subs]))
                    temp = temp.subs(var_subs,Symbol("var('y"+etiqueta.groups()[0]+"',real=True)"))
                f.write("'" + str(Adj(count1)*R*count2) + "'" + ':' + \
                str(temp) + ',')
        f.write( '}')
        f.close()
