# -*- coding: utf-8 -*-
import GMesh as gm
import InOut as IO
import scipy as sp
from scipy import linalg
from scipy import sparse
from scipy.sparse.linalg import spsolve
import os
import sys
import Elements as ele

cwd = os.getcwd()

# -------------------------------------------------------
#     Reading Mesh
# -------------------------------------------------------

#mesh_file = "tube"
mesh_file = "tube-coarse"

global_mesh = gm.GMesh("mesh/" + mesh_file + ".msh")

x_global = global_mesh.X
y_global = global_mesh.Y
ien_global = global_mesh.IEN
nodes_global = len(x_global)
num_ele_global = len(ien_global)

axisym_tri = ele.Linear(x_global, y_global)


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

# -------------------------------------------------------
#     Initial Condition
# -------------------------------------------------------

# Flow
vz = sp.zeros(nodes_global)
vr = sp.zeros(nodes_global)

# Heat
temp_old = sp.zeros(nodes_global)

# -------------------------------------------------------
#     Matrix Assembly
# -------------------------------------------------------

def fem_matrix(_x, _y, _numele, _numnode, _ien):
    kr_global = sp.zeros((_numnode, _numnode), dtype="float64")
    mr_global = sp.zeros((_numnode, _numnode), dtype="float64")
    mr2_global = sp.zeros((_numnode, _numnode), dtype="float64")
    mr_inv_global = sp.zeros((_numnode, _numnode), dtype="float64")
    gxr_global = sp.zeros((_numnode, _numnode), dtype="float64")
    gx = sp.zeros((_numnode, _numnode), dtype="float64")
    gy = sp.zeros((_numnode, _numnode), dtype="float64")
    gyr_global = sp.zeros((_numnode, _numnode), dtype="float64")

    for elem in range(_numele):

        v = [_ien[elem, 0], _ien[elem, 1], _ien[elem, 2]]
        axisym_tri.getAxiSym(v)

        for i_local in range(3):
            i_global = _ien[elem, i_local]
            for j_local in range(3):
                j_global = _ien[elem, j_local]
                kr_global[i_global, j_global] += axisym_tri.kxx[i_local, j_local] \
                                                 + axisym_tri.kyy[i_local, j_local]
                mr_global[i_global, j_global] += axisym_tri.mass[i_local, j_local]
                mr2_global[i_global, j_global] += axisym_tri.mass2[i_local, j_local]
                mr_inv_global[i_global, j_global] += axisym_tri.mass3[i_local, j_local]
                gxr_global[i_global, j_global] += axisym_tri.gx[i_local, j_local]
                gx[i_global, j_global] += axisym_tri.gxm[i_local, j_local]
                gyr_global[i_global, j_global] += axisym_tri.gy[i_local, j_local]
                gy[i_global, j_global] += axisym_tri.gym[i_local, j_local]

    return kr_global, mr_global, mr2_global, mr_inv_global, gxr_global, gyr_global, gx, gy


K, M, M2, M3, Gx, Gy, Gx2, Gy2 = fem_matrix(x_global, y_global, num_ele_global, nodes_global, ien_global)


# ---------------------------------------
# Boundary and Initial Condition
# ---------------------------------------

dirichlet_len = len(global_mesh.dirichlet_points)

vz_a = sp.zeros(nodes_global)
temp_a = sp.zeros(nodes_global)
dpdx = -0.1
for i in range(nodes_global):
    vz_a[i] = -0.25 * dpdx * (1 - y_global[i]**2)
    temp_a[i] = 1


# --------------------------------------
# LHS matrix with Dirichlet BC


# vz = vz_a
# conv = vz * sp.diag(Gx) # + vr * sp.diag(Gy)
# LHS = M * (1./dt) + K * termDiffusivity_fluid + sp.diag(conv)     # + M3 * viscosity_kin + Gy2 + 2*Gx2
LHS = K * termDiffusivity_fluid  # + sp.diag(conv)     # + M3 * viscosity_kin + Gy2 + 2*Gx2
LHS_og = sp.copy(LHS)

cct = sp.zeros(nodes_global)
for i in range(dirichlet_len):
    index = int(global_mesh.dirichlet_points[i][0] - 1)
    value = global_mesh.dirichlet_points[i][1]
    for j in range(nodes_global):
        cct[j] -= value * LHS_og[j, index]
        if j != index:
            LHS[index, j] = 0
            LHS[j, index] = 0
        else:
            LHS[index, j] = 1

# ----------------------------------------------------------
# ---------------------- Loop No Tempo ------------------------


RHS = 1 * cct
for i in range(dirichlet_len):
    index = int(global_mesh.dirichlet_points[i][0]-1)
    RHS[index] = global_mesh.dirichlet_points[i][1]

temp = sp.linalg.solve(LHS, RHS)

# Salvar VTK
vtk_t = IO.InOut(x_global, y_global, ien_global, nodes_global, num_ele_global, temp, None, vz_a
                 , temp_a, temp_a-temp, vz, vr)
vtk_t.saveVTK(cwd+"/results", mesh_file + "-perma")

