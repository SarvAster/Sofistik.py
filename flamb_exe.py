from flamb import *
import numpy as np

cdb_file = r"C:\Users\Côme Delecourt\Desktop\TESTS\Sophistik\barres_exp.cdb"
dat_file = r"C:\Users\Côme Delecourt\Desktop\TESTS\Sophistik\barres_exp.dat"
dll =  r"C:\Program Files\SOFiSTiK\2024\SOFiSTiK 2024\interfaces\64bit\sof_cdb_w-2024.dll"

CDBstatus = CDBinteract(dll)
CDBstatus.open_cdb(cdb_file)
nr, x, y, z = CDBstatus.get_pos()
CDBstatus.close_cdb()
h = max(y)
e = max(x)

V = 5000
H = 2
epsilon = 0.00000001
E = 3.1475870 * 10 ** 7
I = 0.4 * 0.5 **3 / 12
# write get functions for h, E and I
print(h, e)
def analytical_displacement(V, e, h, H, E, I):
    """
    :param V: Force V
    :param e: Eccentricity e
    :param h: Height h
    :param H: Force H
    :param E: Young's Modulus E
    :param I: Moment of Inertia I
    :return: The computed value of w
    """
    k = np.sqrt(V / (E*I))
    term1 = (1 / np.cos(k*h) - 1) * e
    term2 = (np.tan(k*h)/ k - h) * H / V
    w1 = term1 + term2
    w2 = V * e * h**2 / (2 * E * I) + H * h **3 / (3 * E * I)
    return w1, w2

print(analytical_displacement(V, e, h, H, E, I))

iter = Iteration(V, H, epsilon, cdb_file, dat_file, dll)
iter.initialize()
iter.loop()