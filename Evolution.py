from sympy import Symbol, symbols, expand, I, N
from sympy.core.cache import clear_cache
from numpy import conjugate
import shutil
import re as regex

from Metrics import dag, sub_qexpr, op2dict

__all__ = ['init_state', 'sym2num', 'build_ode', 'gsl_main']


def init_state(rho, rho0, basis, dic):
    n = len(basis)
    ini = [0] * (2 * n ** 2)
    tempdic = {}
    for bra in basis:
        bra = dag(bra)
        for ket in basis:
            tempvar = sub_qexpr(bra * rho * ket, dic)
            tempdic[expand((tempvar + conjugate(tempvar)) / 2.)] = expand(
                (bra * rho0 * ket + dag(ket) * rho0 * dag(bra)) / 2.)
            tempdic[expand((tempvar - conjugate(tempvar)) / (2. * I))] = expand(
                (bra * rho0 * ket - dag(ket) * rho0 * dag(bra)) / (2. * I))
    for i in tempdic.keys():
        ini[int(str(i).partition('y')[-1])] = tempdic[i]
    return ini


def sym2num(sol, symbexpr, symbmap=None):
    temp = {}
    L = []
    if symbmap is None:
        symbmap = {s: int(str(s).partition('y')[-1]) for s in symbexpr.atoms(Symbol)}
    for count in range(len(sol)):
        for symbol in symbexpr.atoms(Symbol):
            temp[symbol] = sol[count][symbmap[symbol]]
        L.append(complex(N(symbexpr.subs(temp))))
    return L


def build_ode(rho, rdot, basis, y0=None, t=None, args=(), options={}, filetype=None,
              filename="func.py", dictname="dic.py"):

    if filetype == "GSL":
        D = op2dict(rho, basis)
        f = open(filename, 'w')
        f.write('#include <stdio.h>\n')
        f.write('#include <math.h>\n')
        f.write('#include <complex.h>\n')
        f.write('#include <gsl/gsl_errno.h>\n')
        f.write('#include <gsl/gsl_matrix.h>\n')
        f.write('#include <gsl/gsl_odeiv2.h>\n\n')
        f.write('int func (double t, const double y[], double f[], void *params) {\n')

        W = []
        n = len(basis)
        for k in range(n * n):
            W.append(symbols('Re%d' % k, real=True))
            W.append(symbols('Im%d' % k, real=True))

        L = ['y[' + '%d' % i + ']' for i in range(2 * n * n)]

        lista = {}
        for mi in range(2 * n * n):
            lista[W[mi]] = Symbol(L[mi])

        for i, ele in enumerate(basis):
            bra = dag(ele)
            clear_cache()
            for j, ket in enumerate(basis):
                k = 2 * (i * n + j)
                z = bra * rdot * ket
                zz = sub_qexpr(z, D)
                re = expand((zz + conjugate(zz)) / 2.0, complex=True)
                im = expand((zz - conjugate(zz)) / (2.0 * I), complex=True)
                for var_subs in re.atoms(Symbol):
                    if var_subs in lista:
                        re = re.subs(var_subs, lista[var_subs])
                f.write('    ' + 'f[%d] = ' % k + str(re) + ';\n')
                for var_subs in im.atoms(Symbol):
                    if var_subs in lista:
                        im = im.subs(var_subs, lista[var_subs])
                k = k + 1
                f.write('    ' + 'f[%d] = ' % k + str(im) + ';\n')
        f.write('    ' + 'return GSL_SUCCESS;\n')
        f.write('}\n\n')

        dim = len(L)
        f.write(f'int dim = {str(dim)};\n')

        f = open(dictname, 'w')
        f.write("from sympy import var, I\n")
        f.write('dic=' + '{')
        for bra in basis:
            bra = dag(bra)
            for ket in basis:
                temp = bra * rho * ket
                temp = sub_qexpr(temp, D)
                for var_subs in temp.atoms(Symbol):
                    etiqueta = regex.search(r"\[(.*)\]", str(lista[var_subs]))
                    temp = temp.subs(var_subs, Symbol("var('y" + etiqueta.groups()[0] + "',real=True)"))
                f.write(f"'{str(bra*rho*ket)}': {str(temp)},\n")
        f.write('}')
        f.close()

    elif filetype == "Scipy":
        D = op2dict(rho, basis)
        with open(filename, 'w') as f:
            f.write('from sympy import var, I\n')
            f.write('from numpy import sqrt\n')
            f.write('def f(y, t):\n')

            W = []
            n = len(basis)
            for k in range(n * n):
                W.append(symbols('Re%d' % k, real=True))
                W.append(symbols('Im%d' % k, real=True))

            L = ['y[' + '%d' % i + ']' for i in range(2 * n * n)]

            lista = {}
            for mi in range(2 * n * n):
                lista[W[mi]] = Symbol(L[mi])

            f.write('    ' + 'return [')
            for i, ele in enumerate(basis):
                bra = dag(ele)
                clear_cache()
                for j, ket in enumerate(basis):
                    z = bra * rdot * ket
                    zz = sub_qexpr(z, D)
                    re = expand((zz + conjugate(zz)) / 2.0, complex=True)
                    im = expand((zz - conjugate(zz)) / (2.0 * I), complex=True)
                    for var_subs in re.atoms(Symbol):
                        if var_subs in lista:
                            re = re.subs(var_subs, lista[var_subs])
                    f.write(str(re) + ',\n')
                    for var_subs in im.atoms(Symbol):
                        if var_subs in lista:
                            im = im.subs(var_subs, lista[var_subs])
                    f.write(str(im))
                    if basis.index(ele) * basis.index(ket) == (len(basis) - 1) ** 2:
                        f.write(']\n\n')
                    else:
                        f.write(',\n')

            f.write('dic={')
            for bra in basis:
                bra = dag(bra)
                for ket in basis:
                    temp = bra * rho * ket
                    temp = sub_qexpr(temp, D)
                    for var_subs in temp.atoms(Symbol):
                        etiqueta = regex.search(r"\[(.*)\]", str(lista[var_subs]))
                        temp = temp.subs(var_subs, Symbol("var('y" + etiqueta.groups()[0] + "',real=True)"))
                    f.write("'" + str(bra * rho * ket) + "'" + ':' + str(temp) + ',\n')
            f.write('}')

    elif filetype is None:
        D = op2dict(rho, basis)
        script = 'from scipy.integrate import odeint\n'
        script += 'from sympy import var, I\n'
        script += 'from numpy import sqrt\n'
        script += 'def f(y, t):\n'

        W = []
        n = len(basis)
        for k in range(n * n):
            W.append(symbols('Re%d' % k, real=True))
            W.append(symbols('Im%d' % k, real=True))

        L = ['y[' + '%d' % i + ']' for i in range(2 * n * n)]

        lista = {}
        for mi in range(2 * n * n):
            lista[W[mi]] = Symbol(L[mi])

        script += '    ' + 'return ['
        for i, ele in enumerate(basis):
            bra = dag(ele)
            clear_cache()
            for j, ket in enumerate(basis):
                z = bra * rdot * ket
                zz = sub_qexpr(z, D)
                re = expand((zz + conjugate(zz)) / 2.0, complex=True)
                im = expand((zz - conjugate(zz)) / (2.0 * I), complex=True)
                for var_subs in re.atoms(Symbol):
                    if var_subs in lista:
                        re = re.subs(var_subs, lista[var_subs])
                script += str(re) + ',\n'
                for var_subs in im.atoms(Symbol):
                    if var_subs in lista:
                        im = im.subs(var_subs, lista[var_subs])
                script += str(im)
                if basis.index(ele) * basis.index(ket) == (len(basis) - 1) ** 2:
                    script += ']\n\n'
                else:
                    script += ',\n'

        script += 'dic={'
        for bra in basis:
            bra = dag(bra)
            for ket in basis:
                temp = bra * rho * ket
                temp = sub_qexpr(temp, D)
                for var_subs in temp.atoms(Symbol):
                    etiqueta = regex.search(r"\[(.*)\]", str(lista[var_subs]))
                    temp = temp.subs(var_subs, Symbol("var('y" + etiqueta.groups()[0] + "',real=True)"))
                script += "'" + str(bra * rho * ket) + "'" + ':' + str(temp) + ',\n'
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


def gsl_main(odefile, y0, tspan, step=1000, outfile="funcMain", datafile="dset",
             symbexprs="default", options={}, metadata={}):
    # (sin cambios: no toca DiracObject en absoluto, solo genera texto C)
    h0 = options.get('h0', 1e-6)
    atol = options.get('atol', 1e-6)
    rtol = options.get('rtol', 0.0)
    format = options.get('output_format', "csv")
    solver_method = options.get('solver_method', "rkf45")
    t0, t1 = tspan[0], tspan[1]

    outfile = f"{outfile}.c"
    shutil.copyfile(odefile, outfile)
    if format == "hdf5":
        f = open(outfile, "r")
        contents = f.readlines()
        f.close()
        contents.insert(4, f"#include \"hdf5.h\"\n#define FILENAME \"{datafile}.h5\"\n\n")
        f = open(outfile, "w")
        contents = "".join(contents)
        f.write(contents)
        f.close()
    else:
        pass

    f = open(outfile, 'a')
    f.write('\n')
    f.write('int main(void){\n')
    f.write('    ' + '/*Aqui se colocan los parametros*/ \n')
    f.write('    ' + 'double params = 0; \n')
    f.write('    ' + 'gsl_odeiv2_system sys = {func, NULL, dim, &params};\n')
    f.write('    ' + f'gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_{solver_method}, {str(h0)}, {str(atol)}, {str(rtol)});\n\n')
    f.write('    ' + 'int i;\n')
    f.write('    ' + f'double t = {str(t0)}, t1 = {str(t1)};\n')
    f.write('    ' + 'double y0[] = { ' + str(y0[0]))
    for i in y0[1:]:
        f.write(', ' + str(i))
    f.write(' };\n\n')

    if format == "hdf5":
        if str(symbexprs) == "default":
            f.write('    ' + 'hid_t       file_id, dataset_id, dataspace_id;\n')
            f.write('    ' + 'hsize_t     dims[2];\n')
            f.write('    ' + 'herr_t      status;\n')
            f.write('    ' + 'file_id = H5Fcreate(FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);\n')
            f.write('    ' + f'dims[0] = {str(step)};\n')
            f.write('    ' + 'dims[1] = dim;\n')
            f.write('    ' + 'dataspace_id = H5Screate_simple(2, dims, NULL);\n')
            f.write('    ' + f'dataset_id = H5Dcreate2(file_id, "/{datafile}", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);\n\n')

            if metadata:
                for key, value in metadata.items():
                    f.write('    ' + '{\n')
                    if isinstance(value, list):
                        n_elements = len(value)
                        a_vals = ", ".join([f'"{v}"' for v in value])
                        f.write('    ' + '    ' + f'const char *list_val[{n_elements}] = {{{a_vals}}};\n')
                        f.write('    ' + '    ' + f'hsize_t attr_dims[1] = {{{n_elements}}};\n')
                        f.write('    ' + '    ' + 'hid_t a_space = H5Screate_simple(1, attr_dims, NULL);\n')
                        f.write('    ' + '    ' + 'hid_t type_id = H5Tcopy(H5T_C_S1);\n')
                        f.write('    ' + '    ' + 'H5Tset_size(type_id, H5T_VARIABLE);\n')
                        f.write('    ' + '    ' + f'hid_t attr_id = H5Acreate2(dataset_id, "{key}", type_id, a_space, H5P_DEFAULT, H5P_DEFAULT);\n')
                        f.write('    ' + '    ' + 'H5Awrite(attr_id, type_id, list_val);\n')
                        f.write('    ' + '    ' + 'H5Aclose(attr_id);\n')
                        f.write('    ' + '    ' + 'H5Tclose(type_id);\n')
                        f.write('    ' + '    ' + 'H5Sclose(a_space);\n')
                    else:
                        if isinstance(value, int):
                            h5_type = "H5T_NATIVE_INT"
                            c_type = "int"
                        else:
                            h5_type = "H5T_NATIVE_DOUBLE"
                            c_type = "double"
                        f.write('    ' + '    ' + f'{c_type} val = {value};\n')
                        f.write('    ' + '    ' + 'hid_t a_space = H5Screate(H5S_SCALAR);\n')
                        f.write('    ' + '    ' + f'hid_t attr_id = H5Acreate2(dataset_id, "{key}", {h5_type}, a_space, H5P_DEFAULT, H5P_DEFAULT);\n')
                        f.write('    ' + '    ' + f'H5Awrite(attr_id, {h5_type}, &val);\n')
                        f.write('    ' + '    ' + f'H5Aclose(attr_id);\n')
                        f.write('    ' + '    ' + 'H5Sclose(a_space);\n')
                    f.write('    ' + '}\n')
                f.write('\n')

            f.write('    ' + 'int j, k;\n')
            f.write('    ' + f'double wdata[{str(step)}][dim];\n')
            f.write('    ' + 'for (k = 0; k < dim; k++){\n')
            f.write('    ' + '    ' + 'wdata[0][k] = y0[k];\n')
            f.write('    ' + '}\n')

            f.write('    ' + f'double dt = t1 / {str(step)};\n')
            f.write('    ' + f'for (i = 1; i < {str(step)}; i++)' + '{\n')
            f.write('    ' + '    ' + 'double ti = t + dt;\n')
            f.write('    ' + '    ' + 'int status = gsl_odeiv2_driver_apply (d, &t, ti, y0);\n')
            f.write('    ' + '    ' + 'if (status != GSL_SUCCESS){\n')
            f.write('    ' + '    ' + '    ' + 'printf ("error, return value=%d\\n", status);\n')
            f.write('    ' + '    ' + '    ' + 'break;\n')
            f.write('    ' + '    ' + '}\n')

            f.write('    ' + '    ' + 'for (j = 0; j < dim; j++){\n')
            f.write('    ' + '    ' + '    ' + 'wdata[i][j] = y0[j];\n')
            f.write('    ' + '    ' + '}\n')
        else:
            q7 = len(symbexprs)
            f.write('    ' + 'hid_t       file_id, dataset_id, dataspace_id;\n')
            f.write('    ' + 'hsize_t     dims[2];\n')
            f.write('    ' + 'herr_t      status;\n')
            f.write('    ' + 'file_id = H5Fcreate(FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);\n')
            f.write('    ' + f'dims[0] = {str(step)};\n')
            f.write('    ' + f'dims[1] = {str(q7)};\n')
            f.write('    ' + 'dataspace_id = H5Screate_simple(2, dims, NULL);\n')
            f.write('    ' + f'dataset_id = H5Dcreate2(file_id, "/{datafile}", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);\n\n')

            if metadata:
                for key, value in metadata.items():
                    f.write('    ' + '{\n')
                    if isinstance(value, list):
                        n_elements = len(value)
                        a_vals = ", ".join([f'"{v}"' for v in value])
                        f.write('    ' + '    ' + f'const char *list_val[{n_elements}] = {{{a_vals}}};\n')
                        f.write('    ' + '    ' + f'hsize_t attr_dims[1] = {{{n_elements}}};\n')
                        f.write('    ' + '    ' + 'hid_t a_space = H5Screate_simple(1, attr_dims, NULL);\n')
                        f.write('    ' + '    ' + 'hid_t type_id = H5Tcopy(H5T_C_S1);\n')
                        f.write('    ' + '    ' + 'H5Tset_size(type_id, H5T_VARIABLE);\n')
                        f.write('    ' + '    ' + f'hid_t attr_id = H5Acreate2(dataset_id, "{key}", type_id, a_space, H5P_DEFAULT, H5P_DEFAULT);\n')
                        f.write('    ' + '    ' + 'H5Awrite(attr_id, type_id, list_val);\n')
                        f.write('    ' + '    ' + 'H5Aclose(attr_id);\n')
                        f.write('    ' + '    ' + 'H5Tclose(type_id);\n')
                        f.write('    ' + '    ' + 'H5Sclose(a_space);\n')
                    else:
                        if isinstance(value, int):
                            h5_type = "H5T_NATIVE_INT"
                            c_type = "int"
                        else:
                            h5_type = "H5T_NATIVE_DOUBLE"
                            c_type = "double"
                        f.write('    ' + '    ' + f'{c_type} val = {value};\n')
                        f.write('    ' + '    ' + 'hid_t a_space = H5Screate(H5S_SCALAR);\n')
                        f.write('    ' + '    ' + f'hid_t attr_id = H5Acreate2(dataset_id, "{key}", {h5_type}, a_space, H5P_DEFAULT, H5P_DEFAULT);\n')
                        f.write('    ' + '    ' + f'H5Awrite(attr_id, {h5_type}, &val);\n')
                        f.write('    ' + '    ' + f'H5Aclose(attr_id);\n')
                        f.write('    ' + '    ' + 'H5Sclose(a_space);\n')
                    f.write('    ' + '}\n')
                f.write('\n')

            symbexprs2 = []
            for i in symbexprs:
                expr_str = str(i)
                expr_str = regex.sub(r'y(\d+)', r'y0[\1]', expr_str)
                symbexprs2.append(expr_str)

            f.write('    ' + 'double complex symbexprs2[] = { ' + str(symbexprs2[0]))
            for i in symbexprs2[1:]:
                f.write(', ' + str(i))
            f.write('};\n\n')

            ab = len(symbexprs2)
            f.write('    ' + 'int j, k;\n')
            f.write('    ' + f'double wdata[{str(step)}][{str(ab)}];\n')
            f.write('    ' + f'for (k = 0; k < {str(ab)}; k++)' + '{\n')
            f.write('    ' + '    ' + 'wdata[0][k] = creal(symbexprs2[k]);\n')
            f.write('    ' + '}\n\n')

            f.write('    ' + f'double dt = t1 / {str(step)};\n')
            f.write('    ' + f'for (i = 1; i < {str(step)}; i++)' + '{\n')
            f.write('    ' + '    ' + 'double ti = t + dt;\n')
            f.write('    ' + '    ' + 'int status = gsl_odeiv2_driver_apply (d, &t, ti, y0);\n')
            f.write('    ' + '    ' + 'if (status != GSL_SUCCESS){\n')
            f.write('    ' + '    ' + '    ' + 'printf ("error, return value=%d\\n", status);\n')
            f.write('    ' + '    ' + '    ' + 'break;\n')
            f.write('    ' + '    ' + '}\n')

            f.write('    ' + '    ' + 'double complex symbexprs2[] = { ' + str(symbexprs2[0]))
            for i in symbexprs2[1:]:
                f.write(', ' + str(i))
            f.write('};\n\n')

            f.write('    ' + '    ' + f'int n = {str(ab)};\n\n')
            f.write('    ' + '    ' + 'for (j = 0; j< n; j++){\n')
            f.write('    ' + '    ' + '    ' + 'wdata[i][j] = creal(symbexprs2[j]);\n')
            f.write('    ' + '    ' + '}\n')
        f.write('    ' + '}\n')
        f.write('    ' + 'status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);\n')
        f.write('    ' + 'status = H5Dclose(dataset_id);\n')
        f.write('    ' + 'status = H5Sclose(dataspace_id);\n')
        f.write('    ' + 'status = H5Fclose(file_id);\n')
        f.write('    ' + 'gsl_odeiv2_driver_free (d);\n')
        f.write('    ' + 'return 0;\n')
        f.write('}\n')
        f.close()
    elif format == "csv":
        f.write('    ' + 'FILE *fp;\n')
        f.write('    ' + 'int n;\n')
        f.write('    ' + f'fp = fopen("{datafile}.csv", "a");\n')
        f.write('    ' + 'fprintf (fp, "%.5e", y0[0]);\n')
        f.write('    ' + 'for (n=1 ; n<dim ; n++){\n')
        f.write('    ' + '    ' + 'fprintf (fp, ", %.5e", y0[n]);\n')
        f.write('    ' + '}\n')
        f.write('    ' + 'fclose(fp);\n')
        f.write('    ' + f'double dt = t1 / {str(step)};\n')
        f.write('    ' + f'for (i = 1; i < {str(step)}; i++)' + '{\n')
        f.write('    ' + '    ' + 'double ti = t + dt;\n')
        f.write('    ' + '    ' + 'int status = gsl_odeiv2_driver_apply (d, &t, ti, y0);\n\n')
        f.write('    ' + '    ' + 'if (status != GSL_SUCCESS){\n')
        f.write('    ' + '    ' + '    ' + 'printf ("error, return value=%d\\n", status);\n')
        f.write('    ' + '    ' + '    ' + 'break;\n')
        f.write('    ' + '    ' + '}\n\n')
        f.write('    ' + '    ' + f'fp = fopen("{datafile}.csv", "a");\n')
        f.write('    ' + '    ' + 'int n;\n')
        f.write('    ' + '    ' + 'fprintf (fp, "\\n");\n')
        f.write('    ' + '    ' + 'fprintf (fp, "%.5e", y0[0]);\n')
        f.write('    ' + '    ' + 'for (n=1 ; n<dim ; n++){\n')
        f.write('    ' + '    ' + '    ' + 'fprintf (fp, ", %.5e", y0[n]);\n')
        f.write('    ' + '    ' + '}\n')
        f.write('    ' + '    ' + 'fclose(fp);\n')
        f.write('    ' + '}\n')
        f.write('    ' + 'gsl_odeiv2_driver_free (d);\n')
        f.write('    ' + 'return 0;\n')
        f.write('}\n')
        f.close()