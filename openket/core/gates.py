from diracobject import *
from sympy import I

def X(A = 'default'):
    """
    Pauli-X gate. It works on a single qubit acting like a NOT gate, ie,
    inverts the binary bit changing 0 to 1 and 1 to 0.

    :param A: Operator tag.
    :type A: str, optional
    :return: Returns the Pauli-X gate.
    :rtype: Linear combination of :obj:`Ket`:obj:`Bra` objects

    Example
    ^^^^^^^^^

        .. code-block:: python

            >>> X()
            >>> |0><1| + |1><0|
    
    """
    return Ket(0, A)*Bra(1, A) + Ket(1, A)*Bra(0, A)

def Y(A = 'default'):
    """
    Pauli-Y gate. It works on a single qubit acting as a combination of the Pauli-X and Pauli-Z gates,
    ie, applies both a bit-flip and a phase-flip to the qubit.


    :param A: Operator tag.
    :type A: str, optional
    :return: Returns the Pauli-Y gate.
    :rtype: Linear combination of :obj:`Ket`:obj:`Bra` objects

    Example
    ^^^^^^^^^

        .. code-block:: python

            >>> Y()
            >>> - I|0><1| + I|1><0|
    
    """
    return -I*Ket(0, A)*Bra(1, A) + I*Ket(1, A)*Bra(0, A)

def Z(A = 'default'):
    """
    Pauli-Z gate. It works on a single qubit acting as a phase-flip on it,
    ie, changes the relative phase between the :math:`|0 \\rangle` and :math:`|1 \\rangle` states.


    :param A: Operator tag.
    :type A: str, optional
    :return: Returns the Pauli-Z gate
    :rtype: Linear combination of :obj:`Ket`:obj:`Bra` objects

    Example
    ^^^^^^^^^

        .. code-block:: python

            >>> Z()
            >>> |0><0| - |1><1|

    """
    return Ket(0, A)*Bra(0, A) - Ket(1, A)*Bra(1, A)