# -*- coding: utf-8 -*-
import GMesh as gm
import InOut as IO
import scipy as sp
from scipy import linalg
import matplotlib.pyplot as plt
import Elements as ele
import os

cwd = os.getcwd()

# -------------------------------------------------------
#     Reading Mesh
# -------------------------------------------------------

# mesh_file = "poiseuille"
# mesh_file = "fine"
mesh_file = "poi-var"

savename = mesh_file + "-test"

fluid_mesh = gm.GMesh("mesh/" + mesh_file + "-fld.msh")

x_fluid = fluid_mesh.X
y_fluid = fluid_mesh.Y
ien_fluid = fluid_mesh.IEN
nodes_fluid = len(x_fluid)
num_ele_fluid = len(ien_fluid)

axisym_tri = ele.Linear(x_fluid, y_fluid)

# -------------------------------------------------------
#     Simulation Parameters
# -------------------------------------------------------

# Fluid properties
rho_fluid = 1000
viscosity_din = 1
viscosity_kin = 1

# -------------------------------------------------------
#     Initial Variables
# -------------------------------------------------------

# Flow
psi_last = sp.zeros(nodes_fluid)
vz = sp.zeros(nodes_fluid)
vr = sp.zeros(nodes_fluid)

# -------------------------------------------------------
#     Matrix Assembly
# -------------------------------------------------------


def fem_matrix(_x, _y, _numele, _numnode, _ien):
    k_global = sp.zeros((_numnode, _numnode), dtype="float64")
    mr_global = sp.zeros((_numnode, _numnode), dtype="float64")
    m_global = sp.zeros((_numnode, _numnode), dtype="float64")
    gxr_global = sp.zeros((_numnode, _numnode), dtype="float64")
    gx = sp.zeros((_numnode, _numnode), dtype="float64")
    gy = sp.zeros((_numnode, _numnode), dtype="float64")
    gyr_global = sp.zeros((_numnode, _numnode), dtype="float64")

    for elem in range(_numele):

        v = [_ien[elem, 0], _ien[elem, 1], _ien[elem, 2]]
        axisym_tri.getAxiSym(v)

        ele_radius = (_y[v[0]] + _y[v[1]] + _y[v[2]])/3.

        for i_local in range(3):
            i_global = _ien[elem, i_local]
            for j_local in range(3):
                j_global = _ien[elem, j_local]

                k_global[i_global, j_global]     += ele_radius*axisym_tri.kxx[i_local, j_local]+\
                                                    ele_radius*axisym_tri.kyy[i_local, j_local]+\
                                                    axisym_tri.gy[i_local, j_local]
                mr_global[i_global, j_global]     += ele_radius*axisym_tri.mass[i_local, j_local]
                m_global[i_global, j_global]      += axisym_tri.mass[i_local, j_local]
                gxr_global[i_global, j_global]    += ele_radius*axisym_tri.gx[i_local, j_local]
                gx[i_global, j_global]            += axisym_tri.gx[i_local, j_local]
                gyr_global[i_global, j_global]    += ele_radius*axisym_tri.gy[i_local, j_local]
                gy[i_global, j_global]            += axisym_tri.gy[i_local, j_local]

    return k_global, m_global, mr_global, gxr_global, gyr_global, gx, gy

K, M, Mr, Gxr, Gyr, Gx, Gy = fem_matrix(x_fluid, y_fluid, num_ele_fluid, nodes_fluid, ien_fluid)


MLump = sp.zeros(nodes_fluid)
MinvLump = sp.zeros(nodes_fluid)
for i in range(nodes_fluid):
    for j in range(nodes_fluid):
        MLump[i] += Mr[i, j]
    MinvLump[i] = 1. / MLump[i]

# ---------------------------------------
# Boundary and Initial Condition
# ---------------------------------------

dirichlet_len_fluid = len(fluid_mesh.dirichlet_points)

vz_a = sp.zeros(nodes_fluid)
psi_a = sp.zeros(nodes_fluid)
psi2 = sp.zeros(nodes_fluid)
omega_a = sp.zeros(nodes_fluid)
dpdx = -16
for i in range(nodes_fluid):
    vz_a[i] = -0.25 * dpdx * (1 - y_fluid[i]**2)
    psi_a[i] = -0.25 * dpdx * (0.5 * y_fluid[i]**2 - 0.25 * y_fluid[i]**4)
    omega_a[i] = -0.5 * dpdx * y_fluid[i]


# --------------------------------------
# Psi K matrix with Dirichlet BC -- K_psi

K_psi = K + 2 * Gy
ccpsi = sp.zeros(nodes_fluid)
for i in range(dirichlet_len_fluid):
    index = int(fluid_mesh.dirichlet_points[i][0] - 1)
    value = fluid_mesh.dirichlet_points[i][1]
    for j in range(nodes_fluid):
        ccpsi[j] -= value * K[j, index]
        if j != index:
            K_psi[index, j] = 0
            K_psi[j, index] = 0
        else:
            K_psi[index, j] = 1

F_psi = sp.dot(Mr, omega_a) + ccpsi
for i in range(dirichlet_len_fluid):
    index = int(fluid_mesh.dirichlet_points[i][0] - 1)
    F_psi[index] = fluid_mesh.dirichlet_points[i][1]

psi_last = sp.linalg.solve(K_psi, F_psi)

# Calculo de vz e vr
vz = sp.multiply(MinvLump, sp.dot(Gy, psi_a))
vr = -1.0 * sp.multiply(MinvLump, sp.dot(Gx, psi_a))

# vz = sp.multiply(MinvLump, sp.dot(Gy2, psi_last))
# vr = -1.0 * sp.multiply(MinvLump, sp.dot(Gx2, psi_last))

for i in range(nodes_fluid):
    if y_fluid[i] == max(y_fluid):
        print i
        vz[i] = 0

omega = -1 * sp.multiply(MinvLump,  sp.dot(Gy, vz))


erro = abs(vz_a - vz)
erro_psi = abs(psi_a - psi_last)

# Salvar VTK
vtk_t = IO.InOut(x_fluid, y_fluid, ien_fluid, nodes_fluid, num_ele_fluid, psi_last, omega, vz_a
                 , psi_a, omega_a, vz, vr)
vtk_t.saveVTK(cwd+"/resultsTest", savename + "Test")

