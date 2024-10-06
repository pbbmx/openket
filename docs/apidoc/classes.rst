.. _classes:

***************
Classes
***************


=================================================
OpenKet classes
=================================================


--------------
DiracObject
--------------

.. autoclass:: openket.core.diracobject.DiracObject
    :members:




--------------
Ket
--------------

.. autoclass:: openket.core.diracobject.Ket
    :members:

Examples
^^^^^^^^^
In this example `x` and `v` are two different objects. `x` is a Sympy variable and `v` is the object *Ket*.

    .. code-block:: python

        >>> x = var("x")
        >>> ket_x = Ket(x)
        >>> ket_x
        |x>

In the following example, although `u` and `w` have the same vector tag, they are not treated as the same object. For this kind of tags,
strings as well as Sympy variables may be used.

    .. code-block:: python

        >>> a = Ket(1,"A")
        >>> b = Ket(1,"B")
        >>> a.eig
        1
        >>> b.op
        'B'




--------------
Bra
--------------

.. autoclass:: openket.core.diracobject.Bra
    :members:

Example
^^^^^^^^^
In this example `y` and `z` are two different objects. `y` is a Sympy variable and `z` is the object *Bra*.

    .. code-block:: python
        
        >>> y = var("y")
        >>> bra_y = Bra(y)
        >>> bra_y
        <y|
        >>> bra_y.eig
        y
        >>> bra_y.op
        'default'




--------------
Operator
--------------

.. autoclass:: openket.core.diracobject.Operator
    :members:

Example
^^^^^^^^^
The symbol on the left of the equality does not need to be the same as the argument of Operator, this is just how you are calling it,
where as the argument is how OpenKet tags it.

    .. code-block:: python

        >>> A = Operator("A");
        >>> A
        A
        >>> name = Operator("up")
        >>> name.op
        'up'

It is convenient to have the same tags if no confusion is produced. However the option is always available to name objects the way we want.




------------------------
CreationOperator
------------------------

.. autoclass:: openket.core.diracobject.CreationOperator
    :members:

Example
^^^^^^^^^

    .. code-block:: python

        >>> n = var("n")
        >>> ket_n = Ket(n); bra_n = Bra(n)
        >>> aa = CreationOperator()
        >>> aa*ket_n
        (1 + n)**(1/2)|n+1>
        >>> bra_n*aa
        n**(1/2)<n-1|




------------------------
AnnihilationOperator
------------------------

.. autoclass:: openket.core.diracobject.AnnihilationOperator
    :members:

Example
^^^^^^^^^

    .. code-block:: python

        >>> n = var("n")
        >>> ket_n = Ket(n); bra_n = Bra(n)
        >>> a = AnnihilationOperator()
        >>> a*ket_n
        n**(1/2)|n-1>
        >>> bra_n*a
        (1 + n)**(1/2)<n+1|