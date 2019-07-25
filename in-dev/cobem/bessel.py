# -*- coding: utf-8 -*-
import GMesh as gm
import InOut as IO
import scipy as sp
from scipy import special
import os


cwd = os.getcwd()

# -------------------------------------------------------
#     Reading Mesh
# -------------------------------------------------------

mesh_file = "poiseuille"
savename = "bessel"
mesh = gm.GMesh("mesh/" + mesh_file + "-fld.msh")

x = mesh.X
y = mesh.Y
ien = mesh.IEN
nodes = len(x)
num_ele = len(ien)

temp = sp.zeros(nodes)

# maxsum = 100
# for i in range(nodes):
#     print i, " / ", nodes
#     z = x[i]
#     r = y[i]
#     for n in range(1,maxsum):
#         ln = n*(sp.pi/5.)
#         i0 = sp.special.iv(0, ln)
#         ir = sp.special.iv(0, r*ln)
#         c1 = (2./ln) * (sp.sin(2.5*ln))**2
#         c2 = 2.5 - sp.sin(10*ln)/(4*ln)
#         cn = c1 / (i0 * c2)
#         temp[i] += cn * ir * sp.sin(z*ln)


L = 5
q = 1
k = 1

def gz(_z):
    return (1 - (q*L**2)*(1./(2.*k)) * ((_z/L) - (_z/L)**2)) * sp.sin(ln * _z)

from scipy.integrate import quad

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

