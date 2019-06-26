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

mesh_file = "mesh/poiseuille"

fluid_mesh = gm.GMesh(mesh_file+"-fld.msh")

x_fluid = fluid_mesh.X
y_fluid = fluid_mesh.Y
ien_fluid = fluid_mesh.IEN
nodes_fluid = len(x_fluid)
num_ele_fluid = len(ien_fluid)


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
omega_last = sp.zeros(nodes_fluid)
psi_last = sp.zeros(nodes_fluid)
vz = sp.zeros(nodes_fluid)
vr = sp.zeros(nodes_fluid)

# -------------------------------------------------------
#     Matrix Assembly
# -------------------------------------------------------

def fem_matrix(_x, _y, _numele, _numnode, _ien):
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
    m2_global = sp.zeros((_numnode, _numnode), dtype="float64")
    gx_global = sp.zeros((_numnode, _numnode), dtype="float64")
    gx2 = sp.zeros((_numnode, _numnode), dtype="float64")
    gy_global = sp.zeros((_numnode, _numnode), dtype="float64")

    for elem in range(_numele):
        for i in range(3):
            xx[i] = _x[_ien[elem, i]]
            yy[i] = _y[_ien[elem, i]]

        a[0] = xx[0] * yy[1] - xx[1] * yy[0]        # Change with a[2]?
        a[1] = xx[2] * yy[0] - xx[0] * yy[2]
        a[2] = xx[1] * yy[2] - xx[2] * yy[1]
        area = (a[0] + a[1] + a[2]) / 2.

        radius = (1./3.) * (yy[0] + yy[1] + yy[2])

        b[0] = yy[1] - yy[2]
        b[1] = yy[2] - yy[0]
        b[2] = yy[0] - yy[1]
        c[0] = xx[2] - xx[1]
        c[1] = xx[0] - xx[2]
        c[2] = xx[1] - xx[0]

        for i in range(3):
            for j in range(3):
                k_local[i, j] = (b[i] * b[j] + c[i] * c[j]) * 0.25 * radius * (1. / area)
                gx_local[i, j] = b[j] * (1/6.)
                gy_local[i, j] = c[j] * (1/6.)

        for i_local in range(3):
            i_global = _ien[elem, i_local]
            for j_local in range(3):
                j_global = _ien[elem, j_local]
                k_global[i_global, j_global] += k_local[i_local, j_local]
                m_global[i_global, j_global] += m_local[i_local, j_local] * radius * (area/12.)
                m2_global[i_global, j_global] += m_local[i_local, j_local] * radius * (area/12.)
                gx_global[i_global, j_global] += gx_local[i_local, j_local] * radius
                gx2[i_global, j_global] += gx_local[i_local, j_local] * radius
                gy_global[i_global, j_global] += gy_local[i_local, j_local] * radius

    return k_global, m_global, gx_global, gy_global, m2_global, gx2

K, M, Gx, Gy, M2, gx2 = fem_matrix(x_fluid, y_fluid, num_ele_fluid, nodes_fluid, ien_fluid)







"""

From flowTemp-encit.py:


# --------------------------------------
# Total Matrices (Heat equation)


Alfa_total = sp.zeros(nodes_total)
for i in range(nodes_total):
    if y_total[i] > 0.25:
        Alfa_total[i] = alfa_fld
    else:
        Alfa_total[i] = alfa_solid

K, M, Gx, Gy = fem_matrix(x_total, y_total, num_ele_total, nodes_total, ien_total, Alfa_total)

Mdt = M/dt
# --------------------------------------


# --------------------------------------
# Fluid Region (Psi, Omega)

Alfa_fluid = sp.ones(nodes_fluid)
K_fld, M_fld, Gx_fld, Gy_fld = fem_matrix(x_fluid, y_fluid, num_ele_fluid, nodes_fluid, ien_fluid, Alfa_fluid)

K_fld_ni = K_fld * Ni
Mdt_fld = M_fld/dt

MLump = sp.zeros((nodes_fluid,nodes_fluid))
for i in range(nodes_fluid):
    for j in range(nodes_fluid):
        MLump[i, i] += M_fld[i, j]
MinvLump = linalg.inv(MLump)

K_psi = sp.copy(K_fld)
ccpsi = sp.zeros(nodes_fluid)
# --------------------------------------

# ---------------------------------------
# Condições de contorno e Inicial
# ---------------------------------------

dirichlet_len_total = len(malha_total.dirichlet_points)
dirichlet_len_fluid = len(malha_fluid.dirichlet_points)

Fluid_Boundary = sp.zeros(nodes_fluid)
lista = []
for i in range(nodes_fluid):
    if x_fluid[i] == 0.0 or y_fluid[i] == 0.25 or y_fluid[i] == 0.75 or x_fluid[i] == 2.0:
        Fluid_Boundary[i] = i
    else:
        lista.append(i)
Fluid_Boundary = sp.delete(Fluid_Boundary, lista, axis=0)

num_omega = len(Fluid_Boundary)

# --------------------------------------
# Psi K matrix with Dirichlet BC -- K_psi
for i in range(dirichlet_len_fluid):
    index = int(malha_fluid.dirichlet_points[i][0] - 1)
    value = malha_fluid.dirichlet_points[i][1]
    for j in range(nodes_fluid):
        ccpsi[j] -= value * K_fld[j, index]
        if j != index:
            K_psi[index, j] = 0
            K_psi[j, index] = 0
        else:
            K_psi[index, j] = 1
# --------------------------------------


for i in Fluid_Boundary:
    j = int(i)
    vy[j] = 0.0
    if y_fluid[j] == 0.75:
        vx[j] = 1.0
    if y_fluid[j] == 0.25:
        vx[j] = 0.0

Wz_old = sp.dot(MinvLump, (sp.dot(Gx_fld, vy) - sp.dot(Gy_fld, vx)))

F_psi = sp.dot(M_fld, Wz_old) + ccpsi
for i in range(dirichlet_len_fluid):
    index = int(malha_fluid.dirichlet_points[i][0] - 1)
    F_psi[index] = malha_fluid.dirichlet_points[i][1]

Psi_old = sp.linalg.solve(K_psi, F_psi)


T_a = sp.zeros(nodes_total)
u_a = sp.zeros(nodes_total)
Tan = 20/21.0
for i in range(nodes_total):
    if y_total[i] > 0.25:
        y_linha = y_total[i] - 0.5
        T_a[i] = (Tan/2.) - (Tan * y_linha * 2)
        u_a[i] = 2*y_total[i] - 0.5
    else:
        T_a[i] = ((Tan-1)/0.25)*y_total[i] + 1
        u_a[i] = 0.0


# ----------------------------------------------------------
# ---------------------- Loop No Tempo ------------------------

for t in range(0, tempo):
    print "Solving System " + str((float(t)/(tempo))*100) + "%"

    for i in range(num_omega):
        index = int(Fluid_Boundary[i])
        vy[index] = 0.0
        if y_fluid[index] == 0.75:
            vx[index] = 1.0
        if y_fluid[index] == 0.25:
            vx[index] = 0.0

    for i in range(nodes_total):
        if Vel_points_convert[i] >= 0:
            vx_total[i] = vx[int(Vel_points_convert[i])]
            vy_total[i] = vy[int(Vel_points_convert[i])]
        else:
            vx_total[i] = 0
            vy_total[i] = 0

    # B.C. Vorticidade
    Wcc = sp.dot(MinvLump, (sp.dot(Gx_fld, vy) - sp.dot(Gy_fld, vx)))
    ccomega = sp.zeros(nodes_fluid)

    # Solução de Wz e Psi
    # Conv = sp.dot(sp.diag(vx), Gx) + sp.dot(sp.diag(vy), Gy)
    Conv = vx * Gx_fld + vy * Gy_fld
    LHS_Ni = Mdt_fld + K_fld_ni + Conv
    LHS_omega = sp.copy(LHS_Ni)

    for i in range(num_omega):
        index = int(Fluid_Boundary[i])
        value = Wcc[index]
        for j in range(nodes_fluid):
            ccomega[j] -= value * LHS_Ni[j, index]
            if j != index:
                LHS_omega[index, j] = 0
                LHS_omega[j, index] = 0
            else:
                LHS_omega[index, j] = 1

    F_omega = sp.dot(Mdt_fld, Wz_old) + ccomega  # - sp.dot(Conv, Wz_old)

    for i in range(num_omega):
        index = int(Fluid_Boundary[i])
        F_omega[index] = Wcc[index]

    Wz_new = sp.linalg.solve(LHS_omega, F_omega)

    F_psi = sp.dot(M_fld, Wz_new) + ccpsi
    for i in range(dirichlet_len_fluid):
        index = int(malha_fluid.dirichlet_points[i][0]-1)
        F_psi[index] = malha_fluid.dirichlet_points[i][1]

    Psi_new = sp.linalg.solve(K_psi, F_psi)

    # Temperature OBS: ALMOST DONE, NEED CONVECTIVE TERM
    # Conv_total = sp.dot(sp.diag(vx_total), Gx) + sp.dot(sp.diag(vy_total), Gy)
    Conv_total = vx_total * Gx + vy_total * Gy

    # count=0
    # for i in range(nodes_total):
    #     for j in range(nodes_total):
    #         if Conv_total[i,j] != 0:
    #             count += 1
    #             print Conv_total[i,j], count

    LHS_T = Mdt + K + Conv_total
    LHS_temp = sp.copy(LHS_T)
    cctemp = sp.zeros(nodes_total)
    for i in range(dirichlet_len_total):
        index = int(malha_total.dirichlet_points[i][0]) - 1
        value = malha_total.dirichlet_points[i][1]
        for j in range(nodes_total):
            cctemp[j] -= value * LHS_T[j, index]
            if j != index:
                LHS_temp[index, j] = 0
                LHS_temp[j, index] = 0
            else:
                LHS_temp[index, j] = 1

    F_temp = sp.dot(Mdt, Temp_old) + cctemp
    for i in range(dirichlet_len_total):
        index = int(malha_total.dirichlet_points[i][0]) - 1
        F_temp[index] = malha_total.dirichlet_points[i][1]

    Temp_new = sp.linalg.solve(LHS_temp, F_temp)

    # Salvar VTK
    vtk_t = Io.InOut(x_total, y_total, ien_total, nodes_total, num_ele_total, Temp_old, T_a, u_a, None, None, vx_total, vy_total)
    vtk_t.saveVTK(cwd+"/results", arquivo + str(t+1))

    Psi_old = sp.copy(Psi_new)
    Wz_old = sp.copy(Wz_new)
    Temp_old = sp.copy(Temp_new)

    # Calculo de Vx e Vy
    vx = sp.dot(MinvLump, sp.dot(Gy_fld, Psi_new))
    vy = -1.0 * sp.dot(MinvLump, sp.dot(Gx_fld, Psi_new))

# ----------------- Fim de Loop -------------------
# -------------------------------------------------
"""