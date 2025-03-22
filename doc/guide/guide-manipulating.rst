.. _manipulating:


***********************************
States and operators manipulation
***********************************

In the previous guide section :ref:`basic`, we saw how to use kets, bras and operators,
using the OpenKet classes. Now we will look at some metrics and manipulations on quantum objets.
Also, we will present special basic operators and states, and how they work.


========================
Metrics
========================

------------------------
Conmutator
------------------------

------------------------
Normalize
------------------------

------------------------
Hermitian Conjugate
------------------------

OpenKet is able to calculate the hermitian conjugate of the argument you input, be it a number,
a :obj:`Ket`, a :obj:`Bra`, an :obj:`Operator` or a combination of them, using the :obj:`Adj` function.
For example,

.. testcode:: [basics]
    
    #Definition of bras and kets
    n = var("n")
    ket_n = Ket(n); bra_n = Bra(n)
    #Adj
    print("(I|n>)⁺ =",Adj(ket_n))
    print("(-I<n|)⁺ =",Adj(-sp.I*bra_n))

**Output**:

.. testoutput:: [basics]
    :options: +NORMALIZE_WHITESPACE

    (I|n>)⁺ = <n|
    (-I<n|)⁺ = I|n>

This is particularly useful for computing the density operator of a state, or its norm:

.. testcode:: [basics]
    
    #Definition of the state
    psi = Ket(0) + Ket(1)
    print("We define the state |psi> = ",psi)
    #Adj
    print("|psi>⁺ = <psi| = ",Adj(psi))
    #Norm
    print("<psi|psi> = ",Adj(psi)*psi)
    #Density operator
    print("|psi><psi| = ",psi*Adj(psi))
    print("rho = |psi><psi|/2 = ",psi*Adj(psi)*(1/(Adj(psi)*psi)))

**Output**:

.. testoutput:: [basics]
    :options: +NORMALIZE_WHITESPACE

    We define the state |psi> = |0> + |1>
    |psi>⁺ = <psi| = <0| + <1|
    <psi|psi> = 2
    |psi><psi| = |0><0| + |0><1| + |1><0| + |1><1|
    rho = |psi><psi|/2 = (1/2)|0><0| + (1/2)|0><1| + (1/2)|1><0| + (1/2)|1><1|


------------------------
Trace and Partial Trace
------------------------

OpenKet allows us to trace or partially trace states expressed as a sum of exterior
products. Let's first look into :obj:`TraceOut`. This has 2 arguments: the first one is
the total expression and the second one is the operator tag of the space you will
be tracing out.

.. testcode:: [basics]
    
    #State definition in two spaces
    psi = Ket(0,"H1")*Ket(0,"H2") - Ket(1,"H1")*Ket(1,"H2")
    print("psi = |00> + |11> = ",psi)
    #Proyector operator
    P_psi = psi*Adj(psi)
    print("P_|psi> = |psi><psi| = ",P_psi)
    #TraceOut
    print("Trace of P_|psi> in H1 space: ",TraceOut(P_psi,"H1"))
    print("Trace of P_|psi> in H2 space: ",TraceOut(P_psi,"H2"))

**Output**:

.. testoutput:: [basics]
    :options: +NORMALIZE_WHITESPACE

    psi = |00> + |11> = |0_H1>|0_H2> +  - |1_H1>|1_H2>
    P_|psi> = |psi><psi| = |0_H1>|0_H2><0_H2|<0_H1| +  - |0_H1>|0_H2><1_H2|<1_H1| +  - |1_H1>|1_H2><0_H2|<0_H1| + |1_H1>|1_H2><1_H2|<1_H1|
    Trace of P_|psi> in H1 space: |0_H2><0_H2| + |1_H2><1_H2|
    Trace of P_|psi> in H2 space: |0_H1><0_H1| + |1_H1><1_H1|

The total trace is done with the function :obj:`Trace`. You only need to input the
total expression, `Trace` first finds all the operator tags and will then partially trace
them all out.

.. testcode:: [basics]
    
    #State definition in one space
    eta = Ket(0) - Ket(1)
    print("eta = ",eta)
    P_eta = eta*Adj(eta)
    print("P_|eta> = |eta><eta| = ",P_eta)
    print("Total trace of P_|eta>: ",Trace(P_eta))
    print()
    #State definition in two spaces
    psi = Ket(0,"H1")*Ket(0,"H2") - Ket(1,"H1")*Ket(1,"H2")
    print("psi = |00> + |11> = ",psi)
    P_psi = psi*Adj(psi)
    print("P_|psi> = |psi><psi| = ",P_psi)
    print("Total trace of P_|psi>: ",Trace(P_psi))

**Output**:

.. testoutput:: [basics]
    :options: +NORMALIZE_WHITESPACE

    eta = |0> +  - |1>
    P_|eta> = |eta><eta| = |0><0| +  - |0><1| +  - |1><0| + |1><1|
    Total trace of P_|eta>: 2

    psi = |00> + |11> = |0_H1>|0_H2> +  - |1_H1>|1_H2>
    P_|psi> = |psi><psi| = |0_H1>|0_H2><0_H2|<0_H1| +  - |0_H1>|0_H2><1_H2|<1_H1| +  - |1_H1>|1_H2><0_H2|<0_H1| + |1_H1>|1_H2><1_H2|<1_H1|
    Total trace of P_|psi>: 2


One application of this is to be able to find if we are dealing with a mixed or
a pure state, by computing the square of the density operator and calculating its
total trace.

Suppose we have the pure state :math:`\rho = |\psi \rangle \langle \psi|`, with
:math:`|\psi \rangle = (|0 \rangle - |1 \rangle)/ \sqrt{2}` and we
have the state :math:`\rho_2 = 0.25|0 \rangle \langle 0| + 0.75|1 \rangle \langle 1|`.
First we input this states,

.. testcode:: [basics]
    
    psi = (Ket(0) - Ket(1))*(1/sqrt(2))
    R1 = psi*Adj(psi)
    R2 = 0.25*Ket(0)*Bra(0) + 0.75*Ket(1)*Bra(1)

Both of them have total traces equal to one, but if we square them and calculate their traces,

.. testcode:: [basics]
    
    print("R1 = |psi><psi| = ",R1)
    print("R2 = ",R2)
    print()
    print("Trace(R1*R1) = ",Trace(R1*R1))
    print("Trace(R2*R2) = ",Trace(R2*R2))

**Output**:

.. testoutput:: [basics]
    :options: +NORMALIZE_WHITESPACE

    R1 = |psi><psi| = 0.5|0><0| - 0.5|0><1| - 0.5|1><0| + 0.5|1><1|
    R2 = 0.25|0><0| + 0.75|1><1|

    Trace(R1*R1) = 1
    Trace(R2*R2) = 0.625

This concluding that R2 is not a pure state.


------------------------
Matrix Representation
------------------------

We are able nonetheless to work with the matrix representation of operators. This
is done via :obj:`Qmatrix`. The function has 2 arguments, the first is the expression
(sum of exterior products) and the second is the base in list form.
For example, let's obtain the matrix representation of :math:`\sigma_x` in the computational
basis :math:`\{|0 \rangle, |1 \rangle \}` and :math:`\{|+ \rangle, |- \rangle \}`,
where :math:`|\pm \rangle = (|0 \rangle \pm |1 \rangle)/ \sqrt{2}`.

.. testcode:: [basics]
    
    basis1 =  [Ket(0), Ket(1)]
    basis2 = [(Ket(0)+Ket(1))*(1/sqrt(2)), (Ket(0)-Ket(1))*(1/sqrt(2))]

    print("Qmatrix on basis {|0>,|1>}: ",Qmatrix(X(), base1))
    print("Qmatrix on basis {|+>,|->}: ",Qmatrix(X(), base2))

**Output**:

.. testoutput:: [basics]
    :options: +NORMALIZE_WHITESPACE

    Qmatrix on basis {|0>,|1>}: Matrix([[0, 1], [1, 0]])
    Qmatrix on basis {|+>,|->}: Matrix([[1, 0], [0, -1]])


========================
Built-in operators
========================


========================
Gates
========================
