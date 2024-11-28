.. _basic:


***********************************
Basic operations on Quantum objects
***********************************

========================
Quantum object classes
========================

------------------------
Kets and Bras
------------------------

In quantum mechanics, states are represented as vectors of a Hilbert space. Using Dirac’s notation, these vectors are represented
as :math:`|x \rangle` and are called `kets`. We can define a :obj:`openket.Ket <openket.core.diracobject.Ket>` as follows

.. doctest:: [basics]
    :options: +NORMALIZE_WHITESPACE

    >>> x = var("x")
    >>> xx =  Ket(x)
    >>> print(xx)
    x>

In this example :code:`x` and :code:`xx` are two different objects. :code:`x` is a Sympy variable and :code:`xx` is the object Ket.
The argument given for this object is the vector tag. This tag must be either a number class or a Sympy variable, as in the example.

Vectors corresponding to the dual space are represented by :math:`\langle y|` and are called `bras`. There is also an object
:obj:`openket.Bra <openket.core.diracobject.Bra>` in OpenKet, with exactly the same characteristics as the :obj:`Ket` object,

.. doctest:: [basics]
    :options: +NORMALIZE_WHITESPACE

    >>> y = var("y")
    >>> yy =  Bra(y)
    >>> print(yy)
    <y|


------------------------
Operators
------------------------

An :obj:`openket.Operator <openket.core.diracobject.Operator>` in OpenKet is an object (class) defined as follows

.. doctest:: [basics]
    :options: +NORMALIZE_WHITESPACE

    >>> A = Operator("A")
    >>> print(A)
    A

The symbol on the left of the equality does not need to be the same as the argument of Operator, this is just how you are calling it,
whereas the argument is how OpenKet tags it.


========================
Quantum objects math
========================

------------------------
Eigenvalue equation
------------------------

The eigenvalue equation, :math:`\hat{A} |a \rangle = a |a \rangle` can be computed using OpenKet.

For this, there is another argument (optional) which further specifies the :obj:`Ket` (and :obj:`Bra`). This is the operator tag.
It is necessary to match the names of the operator tags of both the :obj:`Ket` and the :obj:`Operator`

.. doctest:: [basics]
    :options: +NORMALIZE_WHITESPACE

    >>> a = Ket(var("a"), "A")
    >>> A = Operator("A")
    >>> A*a
    a|a_A>


----------------------------------
Basic algebra on bras and kets
----------------------------------

Basic operations may be performed on the :obj:`Kets` and :obj:`Bras` objects, such as:

.. testcode:: [basics]

    b0 = Bra(0); k0 = Ket(0)
    b1 = Bra(1); k1 = Ket(1)
    print('We define: b0 = <0|; k0 = |0>; b1 = <1|; k1 = |1>')

    print('b0*k0 = <0|0> = ',b0 * k0)
    print('k0*b0 = ',k0 * b0)
    print('b0*k1 = <0|1> = ',b0 * k1)
    print('k1*b0 = ',k1 * b0)

**Output**:

.. testoutput:: [basics]
    :options: +NORMALIZE_WHITESPACE

    We define: b0 = <0|; k0 = |0>; b1 = <1|; k1 = |1>
    b0*k0 = <0|0> = 1
    k0*b0 = |0><0|
    b0*k1 = <0|1> = 0
    k1*b0 = |1><0|

-----------------------------
Basic algebra on states
-----------------------------

.. testcode:: [basics]

    #Coefficents definition
    c1 = complex(round(random()*10, 1), round(random()*10, 1))
    c2 = complex(round(random()*10, 1), round(random()*10, 1))
    c3 = complex(round(random()*10, 1), round(random()*10, 1))
    c4 = complex(round(random()*10, 1), round(random()*10, 1))

    #States definition
    psi = c1*k0 + c2*k1
    eta = c3*b0 + c4*b1

    print('psi =', psi)
    print('eta =', eta)
    print('psi*eta = ', psi * eta)
    print('eta*psi = ', eta * psi)

**Output**:

.. testoutput:: [basics]
    :options: +NORMALIZE_WHITESPACE

    psi = (4.9 + 9.6*I)|0> + (2.7 + 4.3*I)|1>
    eta = (2.6 + 1.8*I)<0| + (0.6 + 2.1*I)<1|
    psi*eta =  (2.6 + 1.8*I)*(4.9 + 9.6*I)|0><0| + (0.6 + 2.1*I)*(4.9 + 9.6*I)|0><1| + (2.6 + 1.8*I)*(2.7 + 4.3*I)|1><0| + (0.6 + 2.1*I)*(2.7 + 4.3*I)|1><1|
    eta*psi =  (0.6 + 2.1*I)*(2.7 + 4.3*I) + (2.6 + 1.8*I)*(4.9 + 9.6*I)

