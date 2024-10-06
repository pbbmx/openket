.. _functions:

*********
Functions
*********

=================================================
Metrics on states or operators
=================================================

-------------------
Adj
-------------------

.. autofunction:: openket.core.metrics.Adj

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



-------------------
Normalize
-------------------

.. autofunction:: openket.core.metrics.Normalize



-------------------
Commutator
-------------------

.. autofunction:: openket.core.metrics.Commutator



-------------------
TraceOut
-------------------

.. autofunction:: openket.core.metrics.TraceOut

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



-------------------
Trace
-------------------

.. autofunction:: openket.core.metrics.Trace

Example
^^^^^^^^^
Continuing with the previus example:

    .. code-block:: python

        >>> Trace(R)
        2



=================================================
Visualization
=================================================

-------------------
Qmatrix
-------------------

.. autofunction:: openket.core.metrics.Qmatrix