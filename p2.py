"""
Tests para el port de OpenKet sobre sympy, usando unittest (viene incluido
en Python, no requiere instalar nada aparte de sympy/numpy/scipy que ya
necesita el proyecto).

Cubre: algebra del core, gates de Pauli, metrics (dag/comm/trace/ptrace)
y evolution (generacion de EDOs).

Correr con: python test_openket.py
"""

import unittest
import numpy as np
from sympy import var, sqrt, I, E, cos

from openket_V2 import Ket, Bra, Operator, CreationOperator, AnnihilationOperator
from Metrics import dag, comm, ptrace, trace
from gates import X, Y, Z
from Evolution import build_ode
from Series_v2 import qexp


# ---------------------------------------------------------------------
# Core: creación y aniquilación
# ---------------------------------------------------------------------

class TestCreacionAniquilacion(unittest.TestCase):
    """El comportamiento base viene del docstring original de la libreria,
    asi que si esto se rompe, algo en el motor de evaluacion se rompio."""

    def setUp(self):
        self.n = var("n")
        self.ket_n = Ket(self.n)
        self.bra_n = Bra(self.n)
        self.creador = CreationOperator()
        self.aniquilador = AnnihilationOperator()

    def test_creador_sobre_ket(self):
        self.assertEqual(self.creador * self.ket_n, sqrt(self.n + 1) * Ket(self.n + 1))

    def test_creador_sobre_bra(self):
        self.assertEqual(self.bra_n * self.creador, sqrt(self.n) * Bra(self.n - 1))

    def test_aniquilador_sobre_ket(self):
        self.assertEqual(self.aniquilador * self.ket_n, sqrt(self.n) * Ket(self.n - 1))

    def test_aniquilador_sobre_bra(self):
        self.assertEqual(self.bra_n * self.aniquilador, sqrt(self.n + 1) * Bra(self.n + 1))


class TestCorteMaximo(unittest.TestCase):
    """El parametro `maximum` solo debe frenar la subida del autovalor,
    nunca la bajada -- esto es asimetrico en el diseno original."""

    def test_creador_respeta_el_tope(self):
        creador = CreationOperator("caja", maximum=3)
        self.assertEqual(creador * Ket(2, "caja"), sqrt(3) * Ket(3, "caja"))

    def test_creador_se_anula_en_el_tope(self):
        creador = CreationOperator("caja", maximum=3)
        self.assertEqual(creador * Ket(3, "caja"), 0)

    def test_bajar_el_autovalor_no_respeta_ningun_tope(self):
        # <3|a+ deberia funcionar igual aunque haya un maximum definido,
        # porque bajar el autovalor de un bra nunca esta limitado
        creador = CreationOperator("caja", maximum=3)
        self.assertEqual(Bra(3, "caja") * creador, sqrt(3) * Bra(2, "caja"))


class TestMultiplesEspaciosDeHilbert(unittest.TestCase):

    def test_operador_encuentra_su_ket_saltando_otros_espacios(self):
        """Caso clave para productos tensoriales: un operador que actua
        sobre 'A' debe encontrar su ket aunque haya un ket de otro
        espacio ('B') en medio de la cadena."""
        creador_A = CreationOperator("A")
        resultado = creador_A * Ket(5, "B") * Ket(3, "A")
        self.assertEqual(resultado, 2 * Ket(4, "A") * Ket(5, "B"))

    def test_producto_interno_simple(self):
        self.assertEqual(Bra(3) * Ket(3), 1)
        self.assertEqual(Bra(3) * Ket(4), 0)


# ---------------------------------------------------------------------
# Gates de Pauli
# ---------------------------------------------------------------------

class TestGatesDePauli(unittest.TestCase):

    def test_pauli_x_es_una_compuerta_not(self):
        self.assertEqual(X(), Ket(0) * Bra(1) + Ket(1) * Bra(0))

    def test_pauli_y(self):
        self.assertEqual(Y(), -I * Ket(0) * Bra(1) + I * Ket(1) * Bra(0))

    def test_pauli_z_distingue_signo_entre_0_y_1(self):
        self.assertEqual(Z(), Ket(0) * Bra(0) - Ket(1) * Bra(1))

    def test_pauli_x_al_cuadrado_es_la_identidad_en_el_subespacio(self):
        # X*X deberia devolver al sistema a donde empezo (X es su propia inversa)
        self.assertEqual(X() * X(), Ket(0) * Bra(0) + Ket(1) * Bra(1))


# ---------------------------------------------------------------------
# Metrics: dag, comm, trace, ptrace
# ---------------------------------------------------------------------

class TestDag(unittest.TestCase):

    def test_dag_de_un_numero_complejo(self):
        self.assertEqual(dag(1 + I), 1 - I)

    def test_dag_intercambia_ket_y_bra(self):
        self.assertEqual(dag(Ket(0)), Bra(0))
        self.assertEqual(dag(Bra(1)), Ket(1))

    def test_dag_es_su_propia_inversa_para_operadores(self):
        # Esto fallaba en el codigo original con un NameError (llamaba a
        # una funcion `op()` que no existe en ningun lado del paquete).
        A = Operator("A")
        self.assertEqual(dag(dag(A)), A)


class TestConmutador(unittest.TestCase):

    def test_conmutador_de_dos_operadores_genericos(self):
        A, B = Operator("A"), Operator("B")
        self.assertEqual(comm(A, B), A * B - B * A)


class TestEstadoDeBell(unittest.TestCase):
    """El estado de Bell es el caso de prueba de referencia para todo lo
    que involucre mas de un espacio de Hilbert a la vez."""

    def setUp(self):
        self.bell = (Ket(0, "A") * Ket(0, "B") + Ket(1, "A") * Ket(1, "B")) / sqrt(2)
        self.rho = self.bell * dag(self.bell)

    def test_traza_total_es_uno(self):
        self.assertEqual(trace(self.rho), 1)

    def test_traza_parcial_deja_al_qubit_A_maximamente_mezclado(self):
        rho_A = ptrace(self.rho, "B")
        esperado = (Ket(0, "A") * Bra(0, "A") + Ket(1, "A") * Bra(1, "A")) / 2
        self.assertEqual(rho_A, esperado)


# ---------------------------------------------------------------------
# Evolution: generacion de EDOs
# ---------------------------------------------------------------------

class TestEvolucionTemporal(unittest.TestCase):

    def test_oscilacion_de_rabi_contra_la_solucion_analitica(self):
        """Atomo de dos niveles con H = |0><1| + |1><0|, partiendo de |0>.
        La poblacion en |1> deberia oscilar como sin(t)^2 (Rabi con
        Omega=1). Es el mismo ejemplo que trae el docstring de build_ode."""
        H = Ket(0) * Bra(1) + Ket(1) * Bra(0)
        rho = Operator("R")
        rdot = -I * comm(H, rho)
        basis = [Ket(0), Ket(1)]
        y0 = [1, 0, 0, 0, 0, 0, 0, 0]
        t = np.linspace(0, 10, 200)

        solucion = build_ode(rho=rho, rdot=rdot, basis=basis, y0=y0, t=t, filetype=None)
        poblacion_1 = solucion[:, 6]  # Re<1|R|1>

        error_maximo = np.max(np.abs(poblacion_1 - np.sin(t) ** 2))
        self.assertLess(error_maximo, 1e-5)


# ---------------------------------------------------------------------
# Series de Taylor
# ---------------------------------------------------------------------

class TestSeriesDeTaylor(unittest.TestCase):

    def test_exponencial_de_un_proyector_idempotente(self):
        """e^P para un proyector (P^2=P) tiene forma cerrada conocida:
        1 + (e-1)*P. Con suficientes terminos, deberiamos acercarnos."""
        P = Ket(0) * Bra(0)
        resultado = qexp(P, order=15)
        coeficiente = [t.coef for t in resultado.terms if t.factors == P.factors][0]
        self.assertAlmostEqual(float(coeficiente), float(E) - 1, places=9)

    def test_exponencial_de_un_operador_nilpotente_es_exacta(self):
        # N = |1><0| cumple N^2 = 0 porque <0|1> = 0, asi que la serie se
        # corta sola en el primer termino -- sin aproximacion de por medio.
        N_op = Ket(1) * Bra(0)
        self.assertEqual(N_op ** 2, 0)
        self.assertEqual(qexp(N_op, order=10), 1 + N_op)

    def test_exponencial_imaginaria_de_pauli_x_reproduce_una_rotacion(self):
        # e^(i*theta*X)|0> = cos(theta)|0> + i*sin(theta)|1>, la rotacion
        # de spin de siempre -- y no le dijimos al codigo nada sobre
        # trigonometria, sale sola de la serie de Taylor de sympy mas la
        # idempotencia de X^2=1.
        theta = var("theta")
        estado = qexp(I * theta * X(), order=12) * Ket(0)
        coef_0 = [t.coef for t in estado.terms if t.factors == [Ket(0)]][0]
        esperado = cos(theta).series(theta, 0, 13).removeO()
        self.assertEqual((coef_0 - esperado).expand(), 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)