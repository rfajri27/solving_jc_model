import numpy as np
from qutip import *


def Hamiltonian(omega_a, omega_m, g, n, rwa):
    """
    Mendefinisikan Hamiltonian untuk model Jaynes-Cumming dengan parameter yang sesuai.

    Parameter
    ---------
    omega_a : :int/float
            Merupakan ungkapan frekuensi transisi atom

    omega_m : :int/float
            Merupakan ungkapan frekuensi meden

    g : :int/float
        Merupakan ungkapan parameter kuat interaksi sistem atom-meda

    n : :int
        Merupakan ungkapan parameter jumlah N medan

    rwa : :bool
        Merupakan parameter RWA.
        True : menggunakan RWA
        False : tidak menggunakan RWA

    """

    # Mempersiapkan operator
    a = tensor(destroy(n), qeye(2))  # operator a (anhilation)
    sm = tensor(qeye(n), destroy(2))  # operator sigma-minus

    # Hamiltonian
    if rwa:
        H = omega_m * a.dag() * a + 0.5 * omega_a * commutator(sm.dag(), sm) + g * (a.dag() * sm + a * sm.dag())
    else:
        H = omega_m * a.dag() * a + 0.5 * omega_a * commutator(sm.dag(), sm) + g * (a.dag() + a) * (sm + sm.dag())
    return H
