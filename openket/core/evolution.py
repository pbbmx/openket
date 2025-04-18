from sympy import Symbol, symbols, expand, I, N
from sympy.core.cache import clear_cache
from pylab import conjugate
import re as regex

from .metrics import Adj,Qch,Dictionary

__all__ = ['InitialCondition','SubsSol','Qeq']


def InitialCondition(R, R0, basis, dic):
    """
    This function produces an initial condition in the form of a list (ordered in accordance to the variables :math:`y_i`),
    given the initial condition as a sum of outer products ``R0``.
    The basis ``basis`` has to be given as a list of :obj:`Ket` objects, as well as the dictionary ``dic`` generated by :obj:`openket.Qeq <openket.core.metrics.Qeq>`.
    Also, the density operator ``R`` used in :obj:`openket.Qeq <openket.core.metrics.Qeq>` must be the same.
    
    :param R: Operator involved. Corresponds to the density operator for which we want to calculte the evolution.
    :type R: :obj:`Operator`
    :param R0: Total expression of outer products.
    :type R0: Sum of :obj:`Ket`:obj:`Bra` objects
    :param basis: List of :obj:`Ket` objects containing all the elements of the basis.
    :type basis: list
    :param dic: Dictionary object given by the :obj:`openket.Dictionary <openket.core.metrics.Dictionary>` function.
    :type dic: dict
    :return: Retuns a :math:`2n^2` length list with the initial conditions values, which every two places correspond to a single complex number. :math:`n` is the basis dimension.
    :rtype: list

    Example
    ^^^^^^^^^

    Example with coherent states (quantum optics):
    
        .. code-block:: python
        
            >>> # Define parameters and basis
            >>> n = 5  # Truncation level
            >>> alpha = 1  # Coherent state amplitude
            >>> basis = [Ket(i,"field") for i in range(n)]  # Fock basis
            
            # Create coherent state (initial condition)
            >>> state_alpha = 0
            >>> for i in range(n):
                    state_alpha = state_alpha + ((alpha**i) / math.sqrt(math.factorial(i))) * Ket(i,"field")
            >>> state_alpha = np.exp(-(np.abs(alpha)**2)/2) * state_alpha
            >>> rho0 = state_alpha * Adj(state_alpha)  # Density matrix
            
            # Generate initial conditions for ODE solver input
            >>> init_conditions = InitialCondition(
                    R=Operator("R"), 
                    R0=rho0, 
                    basis=base, 
                    dic=dic
                )
    
            # init_conditions now contains the properly ordered initial values for all math:`y_i` variables
            >>> print(init_conditions)
            [0.367879441171442, 0, 0.367879441171442, 0, 0.260130047511444, 0, 0.150186152955043, 0, 0.0750930764775213, 0, 0.367879441171442, 0, 0.367879441171442, 0, 0.260130047511444, 0, 0.150186152955043, 0, 0.0750930764775213, 0, 0.260130047511444, 0, 0.260130047511444, 0, 0.183939720585721, 0, 0.106197647194831, 0, 0.0530988235974153, 0, 0.150186152955043, 0, 0.150186152955043, 0, 0.106197647194831, 0, 0.0613132401952404, 0, 0.0306566200976202, 0, 0.0750930764775213, 0, 0.0750930764775213, 0, 0.0530988235974153, 0, 0.0306566200976202, 0, 0.0153283100488101, 0]
    
    Note: This example shows how to prepare initial conditions for a coherent state in quantum optics simulations.
    The resulting init_conditions can be used directly with :obj:`scipy.integrate.odeint` or similar ODE solvers.
    """
    n=len(basis)
    ini = [0]*(2*n**2)
    # ini = [0]*(n*(n+1))
    tempdic = {}
    for i in basis:
        for j in basis:
            # if basis.index(i) <= basis.index(j):
                tempvar = Qch(Adj(i)*R*j,dic)
                tempdic[expand( (tempvar + conjugate(tempvar))/2. )] = expand((Adj(i)*R0*j + Adj(j)*R0*i)/2.) #parte real
                tempdic[expand( (tempvar - conjugate(tempvar))/(2.*I) )] = expand((Adj(i)*R0*j - Adj(j)*R0*i)/(2.*I)) #parte imaginaria
    for i in tempdic.keys():
        ini[int(str(i).partition('y')[-1])] = tempdic[i]
    return ini

def SubsSol(sol, var):
    """This function converts symbolic variable expressions into their corresponding time evolution values obtained from an ODE solution.
    It substitutes each symbolic variable in the input expression with its time series data from the ODE solution.
    The function handles both single variables and complex expressions containing multiple variables. For each time point in the solution, it substitutes the 
    variable values and evaluates the expression.

    :param sol: The solution array returned by :obj:`scipy.integrate.odeint`.
            Each column represents the time evolution of a variable (:math:`y0`, :math:`y1`, etc.), with rows corresponding to different time points.
    :type sol: numpy.ndarray
    :param var: symbolic expression of the form :math:`y^k`, given by :obj:`openket.Qch <openket.core.metrics.Qch>`.
    :type var: sympy.expr
    :return: A list containing the numerical evaluated expression values at each time point. The length matches the number of time points in ``sol``.
    :rtype: list
    """
    temp={}; L = []
    if len(var.args) == 0:
        for count1 in range(len(sol)):
            L.append(sol[:, int(str(var).partition('y')[-1])][count1])
            clear_cache()
    else:
        for count1 in range(len(sol)):
            List = var.args
            for count2 in List:
                for count3 in count2.atoms(Symbol):
                    temp[count3] = sol[count1][int(str(count3).partition('y')[-1])]
                    clear_cache()
            L.append(complex(N(var.subs(temp))))
            clear_cache()
    return L

def Qeq(R, Rdot, basis, y0=None, t=None, args=(), options={}, file=None, filename="func", dictname="dic"):
    """
    The main purpose for this function is to study time evolution of an operator, with an expression
    for its time derivative, both given in a specified basis.
    This function creates a file called ``filename`` which defines a function
    in the sintaxis needed to feed `Scipy's ODE solver odeint <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html>`_
    or `GSL's ODE solver gsl_odeiv2.h <https://www.gnu.org/software/gsl/doc/html/ode-initval.html>`_,
    or just applies the :obj:`scipy.integrate.odeint` function with the given parameters.

    :param R: Operator involved. The function turns it in terms of the basis ``basis``.
    :type R: :obj:`Operator`
    :param Rdot: Operator that expresses change over the time, usually the commutator [H,R] where H is the Hamiltonian of the system.
    :type Rdot: :obj:`Operator`
    :param basis: List of :obj:`Ket` objects containing all the elements of the basis.
    :type basis: array
    :param y0: Initial condition on :math:`y`, necessary for the solution only if :code:`file=None`.
    :type y0: array, optional
    :param t: A sequence of time points for which to solve for :math:`y`, necessary if :code:`file=None`.
    :type t: array, optional
    :param args: Extra arguments to pass to :obj:`scipy.integrate.odeint` function (see `documentation <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html>`_).
    :type args: tuple, optional
    :param options:
        * `file_suggestion` (boolean): True if the file created includes a basic code to run the GSL's ODE solver, only available if :code:`file="GSL"`.
        * `h0` (float): The step size to be attempted on the first step, default at 1e-6 for :obj:`gsl_odeiv2.h` and 0.0 for :obj:`scipy.integrate.odeint`.
        * `atol` (float): Absolute tolerance, default at 1e-6 for :obj:`gsl_odeiv2.h` and :code:`None` for :obj:`scipy.integrate.odeint`.
        * `rtol` (float): Relative tolerance, default at 0.0 for :obj:`gsl_odeiv2.h` and :code:`None` for :obj:`scipy.integrate.odeint`.
        * Other Parameters: Other :obj:`scipy.integrate.odeint` parameters (see `documentation <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html>`_).
    :type options: dict, optional
    :param file:
        * :code:`file=None`: the function runs :obj:`scipy.integrate.odeint` with the given parameters.
        * :code:`file="GSL"`: creates a C syntax file to :obj:`gsl_odeiv2.h` solver.
        * :code:`file="Scipy"`: creates a Python syntax file to :obj:`scipy.integrate.odeint` solver.
    :type file: string, optional
    :param filename: The name of the out file with the function.
    :type filename: string, optional
    :param dictname: The name of the dictionay file.
    :type dictname: string, optional
    :return: If :code:`file=None` returns an array with the ODE solution values; other case returns a file with a function written in the sintaxis defined.
            Also returns the dictionary of the base and the variables.
            It first creates a dictionary with the operator ``R`` elements expressed in the basis ``basis``
            using the function :obj:`openket.Dictionary <openket.core.metrics.Dictionary>`, and then it finds the different expressions of time evolution 
            for every matrix element substituting them for variables using :obj:`openket.Qch <openket.core.metrics.Qch>`.
            Finally it writes the function :math:`f` into a file, where the ODE has the form :math:`\\frac{dy}{dt} = f(y,t)`;
            :math:`y` should be understood as a vector, which is defined entry-like as a list, and :math:`t` is the time variable.
            :math:`f` returns also a list, containig the expressions of the derivatives of each entry
            of :math:`y` as entries of another list. 
            :obj:`scipy.integrate.odeint` and :obj:`gsl_odeiv2.h` understand this as a collection of coupled ODE to be solved, specifying
            an initial condition and the time step.
    :rtype: array or file

    Example
    ^^^^^^^^^

    Suppose we want to see the time evolution of the Hamiltonian :math:`\\hat{H} = \\hbar (\\Omega |0 \\rangle \\langle 1| + \\Omega^* |1 \\rangle \\langle 0|)`
    (Two-Level Semiclassical Atom). We assume :math:`\\Omega = 1` for simplicity.
    Following is the script which numerically solves the time evolution of the system described.

        .. code-block:: python
        
            >>> import matplotlib.pyplot as plt
            >>> from numpy import linspace
            >>> from scipy.integrate import odeint

            >>> basis = [Ket(0), Ket(1)]
            >>> rho = Operator("R")
            >>> rho_dot = -I*Commutator(H, rho)
            >>> H = Ket(0)*Bra(1) + Ket(1)*Bra(0)
            
            >>> y0 = [1,0, 0,0, 0,0, 0,0] #[1+i0, 0+i0, 0+i0, 0+i0]
            >>> t = linspace(0,10,1000)

            >>> solution = Qeq(R=rho, Rdot=rho_dot, basis=basis, y0=y0, t=t, file=None)
            >>> plt.plot(t, solution[:,3])
            >>> plt.show()

    Note that the initial condition is given as an array, where the position indicates the bracket term that one wishes
    to give as the initial condition (the dic tells you which is which). As we are working with complex numbers, we represent one number
    with two positions in the list.
    Finally, the command :code:`solution[:,3]` partitions the total solution and only keeps the third entry 
    (:code:`y[3]` corresponding to :math:`\\langle 1| \\rho |1 \\rangle`) in all the lists within `solution`.
    The graph we obtain is the following.

    .. figure:: two_level_system.png
        :scale: 85 %

        Time evolution for :math:`\\langle 1| \\rho |1 \\rangle` with the initial condition :math:`\\rho_{t=0} = |0 \\rangle \\langle 0|`
    
    """
    
    #GSL ODE solver (C language)
    if file == "GSL":  
        #se crea el diccionario para la base b y el operador R
        D = Dictionary(R, basis)
        f = open(filename+'.c', 'w')
        f.write('#include <stdio.h>\n') 
        f.write('#include <complex.h>\n') 
        f.write('#include <gsl/gsl_errno.h>\n')
        f.write('#include <gsl/gsl_matrix.h>\n')
        f.write('#include <gsl/gsl_odeiv2.h>\n\n')
        f.write('int func (double t, const double y[], double f[], void *params) {\n')

        # crea la lista de valores reales e imaginarios W = [Re1,Im1,...,...,Ren*n,Imn*n]
        W = []
        n = len(basis)
        for k in range( n*n ):
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
            # print (i,"/",len(basis))
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
                f.write('    ' + 'f[%d] = ' %k + str(re) + ';\n')
                for var_subs in im.atoms(Symbol):
                    if var_subs in lista:
                        im = im.subs(var_subs,lista[var_subs])
                k = k + 1
                f.write('    ' + 'f[%d] = ' %k + str(im) + ';\n')
        f.write('    ' + 'return GSL_SUCCESS;\n')
        f.write('}\n\n')

        dim = len(L)
        f.write('int main(){\n')
        f.write('    ' + 'int dim = ' + str(dim) + ';\n')
        
        if options.get('file_suggestion', False):
            h0 = options.get('h0', 1e-6)
            atol = options.get('atol', 1e-6)
            rtol = options.get('rtol', 0.0)
            
            lt = len(t) if t is not None else 0
            y0 = str([float(y) for y in y0])[1:-1] if y0 is not None else ""
            t = str(t.tolist())[1:-1] if t is not None else ""

            f.write('    ' + 'gsl_odeiv2_system sys = {func, NULL, dim, NULL}; //{function, jacobian, dimension, params}\n')
            f.write('    ' + f'gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, {str(h0)}, {str(atol)}, {str(rtol)});\n\n')
            f.write('    ' + 'double t[%d] = { '%lt+t+' };\n')
            f.write('    ' + 'double y[%d] = { '%dim+y0+' };\n')
            f.write('    ' + 'for (int i = 0; i <= %d; i++) {\n'%(lt-1))
            f.write('    ' + '    ' + f'double ti = t[i];\n')
            f.write('    ' + '    ' + 'int status = gsl_odeiv2_driver_apply(driver, &t[0], ti, y);\n')
            f.write('    ' + '    ' + 'if (status != GSL_SUCCESS) {\n')
            f.write('    ' + '    ' + '    ' + 'printf("Error en gsl_odeiv2_driver_apply: %s\\n", gsl_strerror(status));\n')
            f.write('    ' + '    ' + '    ' + 'break;\n')
            f.write('    ' + '    ' + '} else {\n')
            f.write('    ' + '    ' + '    ' + 'printf("Step succeeded for t = %.5f\\n", t[0]);\n')
            f.write('    ' + '    ' + '}\n\n')
            f.write('    ' + '    ' + 'printf("t = %.5f\\t", t[0]);\n')
            f.write('    ' + '    ' + 'for (int j = 0; j < dim; j++) {\n')
            f.write('    ' + '    ' + '    ' + 'printf("y[%d] = %.5f\\t", j, y[j]);\n')
            f.write('    ' + '    ' + '}\n')
            f.write('    ' + '    ' + 'printf("\\n");\n')
            f.write('    ' + '}\n')
            f.write('    ' + 'gsl_odeiv2_driver_free(driver);\n')
        f.write('    ' + 'return 0;\n')
        f.write('}')
        f.close()

        # crea el diccionario que relaciona '<i|R|j>' con las y[i]
        f = open(dictname+'.py', 'w')
        f.write('from sympy import var, I\n')
        f.write(dictname + '=' + '{')
        for bra in basis:
            bra = Adj(bra)
            for ket in basis:
                temp = bra*R*ket
                temp = Qch(temp, D)
                for var_subs in temp.atoms(Symbol):
                    etiqueta=regex.search(r"\[(.*)\]",str(lista[var_subs]))
                    temp = temp.subs(var_subs,Symbol("var('y"+etiqueta.groups()[0]+"',real=True)"))
                f.write("'" + str(bra*R*ket) + "'" + ':' + \
                str(temp) + ',\n')
        f.write('}')
        f.close()

    #Scipy's ODE solver (python language)
    elif file == "Scipy":
        D = Dictionary(R, basis)
        with open(filename+".py", 'w') as f:
            f.write('from sympy import var, I\n')
            f.write('from numpy import sqrt\n')
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
                # print (i,"/",len(basis))
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
                        f.write(']\n\n')
                    else:
                        f.write(',\n')

            # crea el diccionario que relaciona '<i|R|j>' con las y[i]
            f.write(dictname + '={')
            for bra in basis:
                bra = Adj(bra)
                for ket in basis:
                    temp = bra*R*ket
                    temp = Qch(temp, D)
                    for var_subs in temp.atoms(Symbol):
                        etiqueta = regex.search(r"\[(.*)\]",str(lista[var_subs]))
                        temp = temp.subs(var_subs,Symbol("var('y"+etiqueta.groups()[0]+"',real=True)"))
                    f.write("'" + str(bra*R*ket) + "'" + ':' + str(temp) + ',\n')
            f.write( '}')

    #Scipy's ODE solver (automatic)
    elif file == None:
        D = Dictionary(R, basis)
        script = 'from scipy.integrate import odeint\n'
        script += 'from sympy import var, I\n'
        script += 'from numpy import sqrt\n'
        script += 'def f(y, t):\n'
        
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
        script += '    ' + 'return ['
        for i,ele in enumerate(basis):
            bra = Adj(ele)
            clear_cache()
            for j,ket in enumerate(basis):
                z = bra*Rdot*ket
                zz = Qch(z, D)
                re = expand((zz+conjugate(zz))/2.0,complex=True)
                im = expand((zz-conjugate(zz))/(2.0*I),complex=True)
                for var_subs in re.atoms(Symbol):
                    if var_subs in lista:
                        re = re.subs(var_subs,lista[var_subs])
                script += str(re)+',\n'
                for var_subs in im.atoms(Symbol):
                    if var_subs in lista:
                        im = im.subs(var_subs,lista[var_subs])
                script += str(im)
                if basis.index(ele)*basis.index(ket)==(len(basis)-1)**2:
                    script += ']\n\n'
                else:
                    script += ',\n'

        # crea el diccionario que relaciona '<i|R|j>' con las y[i]
        script += dictname + '={'
        for bra in basis:
            bra = Adj(bra)
            for ket in basis:
                temp = bra*R*ket
                temp = Qch(temp, D)
                for var_subs in temp.atoms(Symbol):
                    etiqueta = regex.search(r"\[(.*)\]",str(lista[var_subs]))
                    temp = temp.subs(var_subs,Symbol("var('y"+etiqueta.groups()[0]+"',real=True)"))
                script += "'" + str(bra*R*ket) + "'" + ':' + str(temp) + ',\n'
        script += '}\n\n'
        script += f'y0 = {str(y0)}\n'
        script += f't = {t.tolist()}\n'
        script += f"solution = odeint(\
                            func=f,\
                            y0=y0,\
                            t=t,\
                            args={args},\
                            Dfun={options.get('Dfun', None)},\
                            col_deriv={options.get('col_deriv', 0)},\
                            full_output={options.get('full_output', 0)},\
                            ml={options.get('ml', None)},\
                            mu={options.get('mu', None)},\
                            rtol={options.get('rtol', None)},\
                            atol={options.get('atol', None)},\
                            tcrit={options.get('tcrit', None)},\
                            h0={options.get('h0', 0.0)},\
                            hmax={options.get('hmax', 0.0)},\
                            hmin={options.get('hmin', 0.0)},\
                            ixpr={options.get('ixpr', 0)},\
                            mxstep={options.get('mxstep', 0)},\
                            mxhnil={options.get('mxhnil', 0)},\
                            mxordn={options.get('mxordn', 12)},\
                            mxords={options.get('mxords', 5)},\
                            printmessg={options.get('printmessg', 0)},\
                            tfirst={options.get('tfirst', False)}\
                        )"
        context = {}
        exec(script, context)
        return context['solution']
