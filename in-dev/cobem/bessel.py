# -*- coding: utf-8 -*-
import InOut as IO
import scipy as sp
from scipy import special
from scipy.integrate import quad
import matplotlib.pyplot as plt
import os

cwd = os.getcwd()

# -------------------------------------------------------
#  Mesh
# -------------------------------------------------------

L = 5
R = 1
num_x_div = 20
num_y_div = 20

xx = sp.linspace(0, L, num_x_div)
yy = sp.linspace(0, R, num_y_div)
x, y = sp.meshgrid(xx, yy)

nodes = (num_x_div+1) * (num_y_div+1)

# -------------------------------------------------------
#  Parameters
# -------------------------------------------------------

q = 1
k = 1

maxsum = 100

bc = 1


def gz_sin(_z, _ln, _bcFunction = bc):
    return (_bcFunction - (q*L**2)*(1./(2.*k)) * ((_z/L) - (_z/L)**2)) * sp.sin(_ln * _z)


def sin_sqr(_z, _ln):
    return sp.sin(_ln * _z) * sp.sin(_ln * _z)


temperature = sp.zeros(nodes)
for i in range(nodes):
    z = x[i]
    r = y[i]
    for n in range(1, maxsum):
        ln = n*(sp.pi/5.)
        i0 = sp.special.iv(0, ln)
        ir = sp.special.iv(0, r*ln)
        c1 = sp.integrate.quadrature(lambda e: gz_sin(e, ln), 0, L)[0]
        c2 = sp.integrate.quadrature(lambda e: sin_sqr(e, ln), 0, L)[0]
        cn = c1 / (i0 * c2)
        temperature[i] += cn * ir * sp.sin(z*ln)
    temperature[i] += (q*L**2)*(1./(2.*k)) * ((z/L) - (z/L)**2)



plt.pcolor(x, y, temperature)
plt.show()

exit()


def gz(_z):
    return (1 - (q*L**2)*(1./(2.*k)) * ((_z/L) - (_z/L)**2)) * sp.sin(ln * _z)


maxsum = 100
for i in range(nodes):
    print i, " / ", nodes
    z = x[i]
    r = y[i]
    for n in range(1,maxsum):
        ln = n*(sp.pi/5.)
        i0 = sp.special.iv(0, ln)
        ir = sp.special.iv(0, r*ln)
        c1 = quad(gz, 0, L)[0]
        c2 = 2.5 - sp.sin(10*ln)/(4*ln)
        cn = c1 / (i0 * c2)
        temp[i] += cn * ir * sp.sin(z*ln)
    temp[i] += (q*L**2)*(1./(2.*k)) * ((z/L) - (z/L)**2)



# Salvar VTK
vtk_t = IO.InOut(x, y, ien, nodes, num_ele, temp, None, None
                 , None, None, None, None)
vtk_t.saveVTK(cwd+"/results", savename + "Test")

