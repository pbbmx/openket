.. _evolution:


***********************************
Time evolution
***********************************

One important aspect of understanding a certain system is to understand its time
evolution. In quantum mechanics there are many ways of representing how a
system evolves. In this section we will demonstrate the use of OpenKet subroutines
to solve numerically this problem by considering first the Von Neumann equation,

    .. math::

        \dot{\rho} = -\frac{i}{\hbar} [\hat{H}, \rho]

Numerically this is a set of coupled ODE's (ordinary differential equations).
The thing is, we have to find a way to write and name the variables to produce a
script readable by an ODE's solver, and OpenKet helps with this issue.

The problems presented in this user guide are discrete level atoms in special situations.
We face this type of problems by first writing :math:`\hat{H}` as a sum of exterior
products of elements of a basis; then verifying the time evolution of the matrix
elements of :math:`\rho`, namely :math:`\langle i | \rho | j \rangle`, considering
them as variables.

Note: It is convenient to end up with real coefifcients equations, due to the
fact that not all ODE's solvers work with complex numbers, and the point is to make
this as transparent as possible.


========================================
Two-Level Semiclassical Atom Dynamics
========================================

This example models a two-level atom interacting with a classical electromagnetic field
under the rotating wave approximation.

We can describe the atom using a 2 element basis, let this basis be {:math:`|0 \rangle, |1 \rangle`},
where the ground state is :math:`|0 \rangle` and the excited state :math:`|1 \rangle`
Furthermore let the energy of the levels be :math:`\hbar \omega_0` and :math:`\hbar \omega_1`
respectively. Now, the total Hamiltonian may be divided into two parts:
the Hamiltonian of the unperturbed atom :math:`\hat{H_0}`, and the interation part :math:`\hat{V}`.
These can be found to have the following form

    .. math::

        \hat{H_0} = \hbar \omega_0 |0 \rangle \langle 0| + \hbar \omega_1 |1 \rangle \langle 1| \\
        \hat{V} = g E(t) |0 \rangle \langle 1| + g* E(t) |1 \rangle \langle 0|

Where we have put :math:`g = e \langle 0| \vec{r} |1 \rangle \cdot \vec{\epsilon}` and 
:math:`\vec{E}(t) = |\vec{E}(t)|\vec{\epsilon} = E(t) \vec{\epsilon}` with :math:`|\vec{\epsilon}| = 1`.
We assume the field has the form :math:`E(t) = A cos(\nu t)`. It is convenient to change to the interaction picture,
and following the rotating wave approximation (RWA) we end up with the Hamiltonian

    .. math::

        \hat{H} = \hbar (\Omega |0 \rangle \langle 1| + \Omega^* |1 \rangle \langle 0|)

Representing resonant coupling between states, and where :math:`\Omega = \frac{gA}{2 \hbar}`.

--------------------
Code implementation
--------------------

^^^^^^^^^^^^^^^^^^^^
1. System setup
^^^^^^^^^^^^^^^^^^^^

.. testcode:: [basics]

    omega = 1. # Omega = 1 (for simplicity)
    h = 1. # Natural units
    basis = [Ket(0), Ket(1)]  # Two-level basis
    rho = Operator("R")  # Density matrix
    H = h * (omega*Ket(0)*Bra(1) + dag(omega)*Ket(1)*Bra(0))  # Hamiltonian

^^^^^^^^^^^^^^^^^^^^
2. Master Equation
^^^^^^^^^^^^^^^^^^^^

.. testcode:: [basics]

    # Liouville-von Neumann equation (no dissipation)
    rdot = -I*comm(H, rho)

^^^^^^^^^^^^^^^^^^^^^^^^^
3. Initial Conditions
^^^^^^^^^^^^^^^^^^^^^^^^^

.. testcode:: [basics]

    y0 = [1,0, 0,0, 0,0, 0,0]  # Initial state |0⟩⟨0|

^^^^^^^^^^^^^^^^^^^^
4. Time Evolution
^^^^^^^^^^^^^^^^^^^^

In this case, :obj:`openket.build_ode <openket.core.evolution.build_ode>` function provides the ODE's solve using
:obj:`scipy.integrate.odeint` internally, skipping the steps of creating a function and the dictionary
(which in this case do not need to be defined), providing to the function the basis, the initial condition
and the simulation time.

.. testcode:: [basics]

    t = linspace(0,10,1000)  # Simulation time, from 0s to 10s
    solution = build_ode(rho=rho, rdot=rdot, basis=basis, y0=y0, t=t, filetype=None)

^^^^^^^^^^^^^^^^^^^^
5. Visualization
^^^^^^^^^^^^^^^^^^^^

.. testcode:: [basics]

    plt.plot(t, solution[:,3])  # Plot ⟨1|R|1⟩ population
    plt.xlabel('Time (s)')
    plt.ylabel('Excited state population')
    plt.title('Rabi Oscillations in Two-Level System')
    plt.grid(True)
    plt.show()

.. figure:: two_level_system.png
    :scale: 85 %







============================================
Quantum Harmonic Oscillator with Dissipation
============================================

This example demonstrates the dynamics of a quantum harmonic oscillator coupled to a
dissipative environment. We track the time evolution of the expectation values for:
The number of photons :math:`\langle \hat{N} \rangle`, the position :math:`\langle \hat{X} \rangle`
and the momentum :math:`\langle \hat{P} \rangle`.

The Hamiltonian of the system is :math:`\hat{H} = \hbar \omega (\hat{a}^{\dagger} \hat{a} + \frac{1}{2})`,
where :math:`\hat{a}^{\dagger}` and :math:`\hat{a}` are the creation and annihilation operators, and
:math:`\omega` is the oscillator frequency.

The system evolves under the Lindblad master equation:

    .. math::

        \dot{\rho} = -\frac{i}{\hbar} [\hat{H}, \rho] + \frac{\gamma}{2}(2 \hat{a} \rho \hat{a}^{\dagger} - \hat{a}^{\dagger} \hat{a} \rho - \rho \hat{a}^{\dagger} \hat{a})

where :math:`\gamma` is the photon decay rate.

--------------------
Code implementation
--------------------

^^^^^^^^^^^^^^^^^^^^
1. Initialization
^^^^^^^^^^^^^^^^^^^^

.. testcode:: [basics]

    hhbar = 1. # Natural units
    omega = 1. # Oscillator frecuency = 1 (for simplicity)
    gamma = 0.1 # Photon decay rate
    m = 1 # Particle mass = 1 (for simplicity)

    n = 5 # Truncated Fock basis size
    basis = [Ket(i,"field") for i in range(n)] # Basis states |0⟩ to |4⟩

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
2. Operators, Hamiltonian and Master Equation setup
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since OpenKet has already defined :obj:`openket.AnnihilationOperator <openket.core.diracobject.AnnihilationOperator>`
and :obj:`openket.CreationOperator <openket.core.diracobject.CreationOperator>` operators,
we just use them.

.. testcode:: [basics]

    a = AnnihilationOperator("field",n-1)
    aa = CreationOperator("field",n-1)
    rho = Operator("R") # Density matrix
    H = hbar*omega*(aa*a + 1/2) # Hamiltonian

    rdot = -sp.I/hbar * comm(H,rho) + (gamma/2)*(2*a*rho*aa - aa*a*rho - rho*aa*a)

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
3. Initial State: Coherent State
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A coherent state :math:`| \alpha \rangle = e^{-\frac{|\alpha|^2}{2}} \sum_{n=0}^{\infty} \frac{\alpha^n}{\sqrt{n!}} |n \rangle`,
with :math:`\alpha = 1` (mean photon number :math:`\langle \hat{N} \rangle (0)=1`)

.. testcode:: [basics]

    alpha = 1
    state_alpha = 0
    for i in range(n):
            state_alpha = state_alpha + ((alpha**2) / math.sqrt(math.factorial(i))) * Ket(i,"field")
    state_alpha = np.exp(-(np.abs(alpha)**2)/2) * state_alpha
    rho0 = state_alpha * dag(state_alpha) # Density matrix ρ(0) = |α⟩⟨α|

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
4. Time Evolution (Scipy's ODE solver)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :obj:`openket.build_ode <openket.core.evolution.build_ode>` function provides the coupled ODE's that must
be solved to find the solution, solving the issue mentioned at the beginning of the chapter, not forgetting that since
we want to use Scipy's ODE solver, we must specify :code:`filetype="Scipy"`.
Also, we are using the :obj:`openket.init_state <openket.core.evolution.init_state>`
function to represent the coherent states, used as a initial conditions, as a list of values.
Both are necessary to be fed into :obj:`scipy.integrate.odeint`, or any other ODE's solver.

.. testcode:: [basics]

    build_ode(rho=rho, rdot=rdot, basis=basis, filetype="Scipy", filename="func.py")
    from func import dic,f
    init_conditions = init_state(rho=rho, rho0=rho0, basis=basis, dic=dic)
    t = np.linspace(0,50,500) # Simulation time, from 0s to 50s
    solution = odeint(f, init_conditions, t)

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
5. Expectation Values (Scipy's ODE solver)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To compute :math:`\langle \hat{N} \rangle`, :math:`\langle \hat{X} \rangle` and
:math:`\langle \hat{P} \rangle`, we must define the :math:`\hat{N}`, :math:`\hat{X}`
and :math:`\hat{P}` operators, as:

* :math:`\hat{N} = \hat{a}^{\dagger} \hat{a}`
* :math:`\hat{X} = \sqrt{\frac{\hbar}{2 m \omega}} (\hat{a} + \hat{a}^{\dagger})`
* :math:`\hat{P} = i \sqrt{\frac{\hbar m \omega}{2}} (\hat{a}^{\dagger} - \hat{a})`

Then, using the property :math:`\langle \hat{A} \rangle = Tr[\rho \hat{A}]`,
we can find the expression of the trace using :obj:`openket.trace <openket.core.metrics.trace>`
in terms of the given basis, and then replace the :math:`\langle i | R | j \rangle`
objects with their corresponding variables :math:`y^k`, using :obj:`openket.sub_qexpr <openket.core.metrics.sub_qexpr>`.
Now, with :obj:`openket.sym2num <openket.core.evolution.sym2num>`,
we can find the correct numerical value in the density matrix with the correct position in the list.
As each row of :code:`solution` corresponds to density matrix :math:`\rho` at a given time,
we get a list with the values at each time :math:`t`.

.. testcode:: [basics]

    N = aa * a
    x = np.sqrt(hbar/(2*m*omega)) * (a + aa)
    p = sp.I * np.sqrt(hbar*m*omega/2) * (aa - a)

    N_symb = sub_qexpr(qexpr=trace(rho * N, basis=basis), dic=dic)
    x_symb = sub_qexpr(qexpr=trace(rho * x, basis=basis), dic=dic)
    p_symb = sub_qexpr(qexpr=trace(rho * p, basis=basis), dic=dic)

    N_expect = sym2num(sol=solution, symbexpr=N_symb)
    x_expect = sym2num(sol=solution, symbexpr=x_symb)
    p_expect = sym2num(sol=solution, symbexpr=p_symb)


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
6. Time Evolution and Expectation Values (GSL's ODE solver)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Again, the :obj:`openket.build_ode <openket.core.evolution.build_ode>` function provides the coupled ODE's that must
be solved to find the solution, but this time in C code written in another file.
As before, we are using the :obj:`openket.init_state <openket.core.evolution.init_state>`
function to represent the coherent states (initial conditions) as a list of values, but now
we must set :code:`filetype="GSL"`.

.. testcode:: [basics]

    filename = "gslode.c"
    build_ode(rho=rho, rdot=rdot, basis=base, filetype="GSL", filename=filename)
    from dic import dic
    init_conditions = init_state(rho=rho, rho0=rho0, basis=base, dic=dic)
    t0 = 0, t1 = 50, step = 500
    t = np.linspace(t0,t1,step)


We define the :math:`\hat{N}`, :math:`\hat{X}` and :math:`\hat{P}` operators and now we have all the elements
to feed :obj:`openket.gsl_main <openket.core.evolution.gsl_main>` function, obtaining a C file named `"gslode_main.c"`.

.. testcode:: [basics]

    N = aa * a
    x = np.sqrt(hbar/(2*m*omega)) * (a + aa)
    p = I * np.sqrt(hbar*m*omega/2) * (aa - a)

    N_symb = sub_qexpr(qexpr=trace(rho * N, basis=base), dic=dic)
    x_symb = sub_qexpr(qexpr=trace(rho * x, basis=base), dic=dic)
    p_symb = sub_qexpr(qexpr=trace(rho * p, basis=base), dic=dic)

    gsl_main(odefile=filename,
        outfile="gslode_main.c",
        y0=init_conditions,
        tspan=(t0,t1),
        step=step,
        symbexprs=[N_symb,x_symb,p_symb],
        options={ "output_format": "hdf5" })

To compile the C code, we simply type in the terminal

.. code-block:: bash

    gcc gslode_main.c -lgsl -lgslcblas -lhdf5 -lm
    ./a.out


And we obtain a dataset with HDF5 format (in this case) that we can read with the following lines:

.. testcode:: [basics]

    with h5py.File('dset.h5', 'r') as f:
        N_expect = f['dset'][:,0]
        x_expect = f['dset'][:,1]
        p_expect = f['dset'][:,2]

^^^^^^^^^^^^^^^
7. Plotting
^^^^^^^^^^^^^^^

.. testcode:: [basics]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Subplot 1: <N>
    ax1.plot(t, N_expect, 'r-', label='<N>')
    ax1.set_ylabel('Photon number')
    ax1.legend()
    ax1.grid(True)

    # Subplot 2: <X> and <P>
    ax2.plot(t, x_expect, 'b-', label='<X>')
    ax2.plot(t, p_expect, 'g--', label='<P>')
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Expectation value')
    ax2.legend()
    ax2.grid(True)

    plt.suptitle('Dynamics of a Dissipative Quantum Harmonic Oscillator')
    plt.tight_layout()
    plt.show()

.. figure:: dissipative_QHO.png
    :scale: 85 %