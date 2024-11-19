.. _basic:

******************
Introduction
******************

This is an introductory manual to the Python library OpenKet. This library is a tool for manipulating objects from quantum mechanics,
such as vectors in the form of kets, operators, etc.
It is free software under GNU-GPL v3 license, and uses other free software libraries, such as Sympy, Numpy and Scipy.
We also provide the API documentation in :ref:`apidoc`.



.. warning:: Do not run OpenKet from the installation directory.

Before typing the examples presented in the following sections, it is necessary to import OpenKet

.. code-block:: Python

   from openket import *

This will load all of the user available functions. Often, we also need to import the NumPy, Simpy and SciPy libraries with:

.. code-block:: Python

    import numpy as np
    import simpy as sp
    import scipy as sc

In the rest of the documentation, functions are written using `qutip.module.function()` notation which links to the corresponding
function in the OpenKet API: :ref:`apidoc`. However, in calling `import *`, we have already loaded all of the OpenKet modules.
Therefore, we will only need the function name and not the complete path when calling the function from the interpreter prompt,
Python script, or Jupyter notebook.



========================
Kets and Bras
========================

In quantum mechanics, states are represented as vectors of a Hilbert space. Using Dirac’s notation, these vectors are represented
as :math:`|x \rangle` and are called kets. OpenKet understands kets as an object (class) with specific properties or methods. We can
define a ket as follows

.. testcode:: [basics]

    x = var("x")
    xx =  Ket(x)
    print(xx)

**Output**:

.. testoutput:: [basics]
    :options: +NORMALIZE_WHITESPACE

    |x>