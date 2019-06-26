# -*- coding: utf-8 -*-
import GMesh as gm
import InOut as IO
import scipy as sp
from scipy import linalg
from scipy import sparse
from scipy.sparse.linalg import spsolve
import os
import sys

cwd = os.getcwd()

# -------------------------------------------------------
#     Reading Mesh
# -------------------------------------------------------

mesh_file = "mesh/test"

global_mesh = gm.GMesh(mesh_file+".msh")

x_global = global_mesh.X
y_global = global_mesh.Y
ien_global = global_mesh.IEN
nodes_global = len(x_global)
num_ele_global = len(ien_global)

fluid_mesh = gm.GMesh(mesh_file+"-fld.msh")

x_fluid = fluid_mesh.X
y_fluid = fluid_mesh.Y
ien_fluid = fluid_mesh.IEN
nodes_fluid = len(x_fluid)
num_ele_fluid = len(ien_fluid)

Vel_points_convert = sp.zeros(nodes_global) -1
for i in range(nodes_global):
    for j in range(nodes_fluid):
        if (x_global[i] == x_fluid[j]) and (y_global[i] == y_fluid[j]):
            Vel_points_convert[i] = j

# -------------------------------------------------------
#     Simulation Parameters
# -------------------------------------------------------

dt = 0.01
dt_inv = 1. / dt
time = 1000

# Fluid properties
rho_fluid = 1000
viscostity_din = 0.8e-03
viscostity_kin = viscostity_din / rho_fluid
termCondutivity_fluid = 0.6089
spHeat_fluid = 4137.9
termDiffusivity_fluid = termCondutivity_fluid / (rho_fluid * spHeat_fluid)

# Solid properties
rho_solid = 2000
termCondutivity_solid = 10
spHeat_solid = 10
termDiffusivity_solid = termCondutivity_solid / (rho_solid * spHeat_solid)

# -------------------------------------------------------
#     Initial Condition
# -------------------------------------------------------

# Flow
omega_last = sp.zeros(nodes_fluid)
psi_last = sp.zeros(nodes_fluid)
vz = sp.zeros(nodes_fluid)
vr = sp.zeros(nodes_fluid)

# Heat
temp_last = sp.zeros(nodes_global)

# -------------------------------------------------------
#     Matrix Assembly
# -------------------------------------------------------

def fem_matrix(_x, _y, _numele, _numnode, _ien, _alfa):
    k_local = sp.zeros((3, 3), dtype="float64")
    m_local = sp.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype="float64")
    gx_local = sp.zeros((3, 3), dtype="float64")
    gy_local = sp.zeros((3, 3), dtype="float64")
    a = sp.zeros(3, dtype="float64")
    b = sp.zeros(3, dtype="float64")
    c = sp.zeros(3, dtype="float64")
    yy = sp.zeros(3, dtype="float64")
    xx = sp.zeros(3, dtype="float64")
    k_global = sp.zeros((_numnode, _numnode), dtype="float64")
    m_global = sp.zeros((_numnode, _numnode), dtype="float64")
    gx_global = sp.zeros((_numnode, _numnode), dtype="float64")
    gy_global = sp.zeros((_numnode, _numnode), dtype="float64")

    for elem in range(_numele):
        alfa = 0
        for i in range(3):
            xx[i] = _x[_ien[elem, i]]
            yy[i] = _y[_ien[elem, i]]
            alfa += _alfa[_ien[elem, i]]

        alfa = alfa * (1/3.0)
        a[0] = xx[0] * yy[1] - xx[1] * yy[0]
        a[1] = xx[2] * yy[0] - xx[0] * yy[2]
        a[2] = xx[1] * yy[2] - xx[2] * yy[1]
        Area = (a[0] + a[1] + a[2]) / 2.

        b[0] = yy[1] - yy[2]
        b[1] = yy[2] - yy[0]
        b[2] = yy[0] - yy[1]
        c[0] = xx[2] - xx[1]
        c[1] = xx[0] - xx[2]
        c[2] = xx[1] - xx[0]

        for i in range(3):
            for j in range(3):
                k_local[i, j] = (b[i] * b[j] + c[i] * c[j]) / (4 * Area)
                gx_local[i,j] = b[j] * (1/6.)
                gy_local[i,j] = c[j] * (1/6.)

        for i_local in range(3):
            i_global = _ien[elem, i_local]
            for j_local in range(3):
                j_global = _ien[elem, j_local]
                k_global[i_global, j_global] += k_local[i_local, j_local]*alfa
                m_global[i_global, j_global] += m_local[i_local, j_local]* (Area/12.)
                gx_global[i_global, j_global] += gx_local[i_local, j_local]
                gy_global[i_global, j_global] += gy_local[i_local, j_local]

    return  k_global, m_global, gx_global, gy_global