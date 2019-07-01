# -*- coding: utf-8 -*-
import GMesh as gm
import InOut as IO
import scipy as sp
from scipy import linalg
import os

cwd = os.getcwd()

# -------------------------------------------------------
#     Reading Mesh
# -------------------------------------------------------

mesh_file = "poiseuille"
# mesh_file = "fine"

fluid_mesh = gm.GMesh("mesh/" + mesh_file + "-fld.msh")

x_fluid = fluid_mesh.X
y_fluid = fluid_mesh.Y
ien_fluid = fluid_mesh.IEN
nodes_fluid = len(x_fluid)
num_ele_fluid = len(ien_fluid)


# -------------------------------------------------------
#     Simulation Parameters
# -------------------------------------------------------

dt = 0.005
time = 400

# Fluid properties
rho_fluid = 1000
viscosity_din = 1
viscosity_kin = 1

# -------------------------------------------------------
#     Initial Variables
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
    m3_global = sp.zeros((_numnode, _numnode), dtype="float64")
    gx_global = sp.zeros((_numnode, _numnode), dtype="float64")
    gx2 = sp.zeros((_numnode, _numnode), dtype="float64")
    gy2 = sp.zeros((_numnode, _numnode), dtype="float64")
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
                k_local[i, j] = (b[i] * b[j] + c[i] * c[j]) * 0.25 * radius * (1./ area)
                gx_local[i, j] = b[j] * (1/6.)
                gy_local[i, j] = c[j] * (1/6.)

        for i_local in range(3):
            i_global = _ien[elem, i_local]
            for j_local in range(3):
                j_global = _ien[elem, j_local]
                k_global[i_global, j_global] += k_local[i_local, j_local]
                m_global[i_global, j_global] += m_local[i_local, j_local] * radius * (area/12.)
                m2_global[i_global, j_global] += m_local[i_local, j_local] * radius**2 * (area/12.)
                m3_global[i_global, j_global] += m_local[i_local, j_local] * (1./radius**2) * (area/12.)
                gx_global[i_global, j_global] += gx_local[i_local, j_local] * radius
                gx2[i_global, j_global] += gx_local[i_local, j_local]
                gy_global[i_global, j_global] += gy_local[i_local, j_local] * radius
                gy2[i_global, j_global] += gy_local[i_local, j_local]

    return k_global, m_global, gx_global, gy_global, m2_global, gx2, gy2, m3_global


K, M, Gx, Gy, M2, Gx2, Gy2, M3 = fem_matrix(x_fluid, y_fluid, num_ele_fluid, nodes_fluid, ien_fluid)

Mdt = M/dt
K_ni = K * viscosity_kin


MLump = sp.zeros(nodes_fluid)
MLump3 = sp.zeros(nodes_fluid)
MinvLump = sp.zeros(nodes_fluid)
for i in range(nodes_fluid):
    for j in range(nodes_fluid):
        MLump[i] += M[i, j]
        MLump3[i] += M3[i, j]
    MinvLump[i] = 1. / MLump[i]



# ---------------------------------------
# Boundary and Initial Condition
# ---------------------------------------

dirichlet_len_fluid = len(fluid_mesh.dirichlet_points)

Fluid_Boundary_in = sp.zeros(nodes_fluid)
Fluid_Boundary_out = sp.zeros(nodes_fluid)
Fluid_Boundary_wall = sp.zeros(nodes_fluid)
Fluid_Boundary_axis = sp.zeros(nodes_fluid)

Fluid_Boundary = sp.zeros(nodes_fluid)

list1 = []
list2 = []
list3 = []
list4 = []
list5 = []

for i in range(nodes_fluid):
    if x_fluid[i] == 0:
        Fluid_Boundary_in[i] = i
    else:
        list1.append(i)

    if x_fluid[i] == 5:
        Fluid_Boundary_out[i] = i
    else:
        list2.append(i)

    if y_fluid[i] == 1:
        Fluid_Boundary_wall[i] = i
    else:
        list3.append(i)

    if y_fluid[i] == 0:
        Fluid_Boundary_axis[i] = i
    else:
        list4.append(i)

    if x_fluid[i] == 0 or x_fluid[i] == 5 \
            or y_fluid[i] == 0 or y_fluid[i] == 1:
        Fluid_Boundary[i] = i
    else:
        list5.append(i)


Fluid_Boundary_in = sp.delete(Fluid_Boundary_in, list1, axis=0)
Fluid_Boundary_out = sp.delete(Fluid_Boundary_out, list2, axis=0)
Fluid_Boundary_wall = sp.delete(Fluid_Boundary_wall, list3, axis=0)
Fluid_Boundary_axis = sp.delete(Fluid_Boundary_axis, list4, axis=0)

Fluid_Boundary = sp.delete(Fluid_Boundary, list5, axis=0)


num_omega = len(Fluid_Boundary)


vz_a = sp.zeros(nodes_fluid)
psi_a = sp.zeros(nodes_fluid)
omega_a = sp.zeros(nodes_fluid)
dpdx = -16
for i in range(nodes_fluid):
    vz_a[i] = -0.25 * dpdx * (1 - y_fluid[i]**2)
    psi_a[i] = -0.25 * dpdx * (0.5 * y_fluid[i]**2 - 0.25 * y_fluid[i]**4)
    omega_a[i] = -0.5 * dpdx * y_fluid[i]

# --------------------------------------
# Psi K matrix with Dirichlet BC -- K_psi

K_psi = 1*K + 2 * Gy2
# K_psi = 1*K + 2 * Gy
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
# --------------------------------------

for i in Fluid_Boundary:
    j = int(i)
    vr[j] = 0.0
    if x_fluid[j] == 0:
        # vz[j] = 1.0
        vz[j] = vz_a[j]
    if y_fluid[j] == 1:
        vz[j] = 0
        vr[j] = 0
   # if y_fluid[j] == 0:
    #    vz[j] = 4

omega_last = sp.multiply(MinvLump, (sp.dot(Gx, vr) - sp.dot(Gy, vz)))

F_psi = sp.dot(M2, omega_last) + ccpsi
# F_psi = sp.dot(M, omega_last) + ccpsi
for i in range(dirichlet_len_fluid):
    index = int(fluid_mesh.dirichlet_points[i][0] - 1)
    F_psi[index] = fluid_mesh.dirichlet_points[i][1]

psi_last = sp.linalg.solve(K_psi, F_psi)


# ----------------------------------------------------------
# ---------------------- Loop No Tempo ------------------------

for t in range(0, time):
    print "Solving System " + str((float(t)/time)*100) + "%"

    # vr = sp.zeros(nodes_fluid)

    for i in range(num_omega):
        index = int(Fluid_Boundary[i])
        vr[index] = 0
        if y_fluid[index] == 1:
            vz[index] = 0.
            vr[index] = 0.
        if x_fluid[index] == 0:
            # vz[index] = 1.
            vz[index] = vz_a[index]
        if y_fluid[index] == 0:
            vz[index] = 4.

    # B.C. Vorticidade
    W_in = sp.multiply(MinvLump, (sp.dot(Gx, vr) - sp.dot(Gy, vz)))
    W_axis = 0.
    W_wall = W_in
    ccomega = sp.zeros(nodes_fluid)

    # Solução de Wz e Psi
    Conv = vz * sp.diag(Gx) + vr * sp.diag(Gy)

    LHS_Ni = Mdt + K_ni + sp.diag(Conv) + MLump3*viscosity_kin
    LHS_omega = sp.copy(LHS_Ni)

    for i in range(len(Fluid_Boundary_in)):
       index = int(Fluid_Boundary_in[i])
       value = 8*y_fluid[index]
       for j in range(nodes_fluid):
           ccomega[j] -= value * LHS_Ni[j, index]
           if j != index:
               LHS_omega[index, j] = 0
               LHS_omega[j, index] = 0
           else:
               LHS_omega[index, j] = 1

    for i in range(len(Fluid_Boundary_axis)):
        index = int(Fluid_Boundary_axis[i])
        value = W_axis
        for j in range(nodes_fluid):
            ccomega[j] -= value * LHS_Ni[j, index]
            if j != index:
                LHS_omega[index, j] = 0
                LHS_omega[j, index] = 0
            else:
                LHS_omega[index, j] = 1

    for i in range(len(Fluid_Boundary_wall)):
        index = int(Fluid_Boundary_wall[i])
        value = 8
        for j in range(nodes_fluid):
            ccomega[j] -= value * LHS_Ni[j, index]
            if j != index:
                LHS_omega[index, j] = 0
                LHS_omega[j, index] = 0
            else:
                LHS_omega[index, j] = 1

    F_omega = sp.dot(Mdt, omega_last) + ccomega  # - sp.dot(Conv, omega_last)

    for i in range(len(Fluid_Boundary_in)):
       index = int(Fluid_Boundary_in[i])
       F_omega[index] = 8*y_fluid[index]

    for i in range(len(Fluid_Boundary_wall)):
        index = int(Fluid_Boundary_wall[i])
        F_omega[index] = 8

    for i in range(len(Fluid_Boundary_axis)):
        index = int(Fluid_Boundary_axis[i])
        F_omega[index] = W_axis


    omega = sp.linalg.solve(LHS_omega, F_omega)

    F_psi = sp.dot(M2, omega) + ccpsi
    # F_psi = sp.dot(M, omega) + ccpsi
    for i in range(dirichlet_len_fluid):
        index = int(fluid_mesh.dirichlet_points[i][0]-1)
        F_psi[index] = fluid_mesh.dirichlet_points[i][1]

    psi = sp.linalg.solve(K_psi, F_psi)


    # Salvar VTK
    vtk_t = IO.InOut(x_fluid, y_fluid, ien_fluid, nodes_fluid, num_ele_fluid, psi, omega, vz_a
                     , psi_a, omega_a, vz, vr)
    vtk_t.saveVTK(cwd+"/results", mesh_file + str(t+1))

    Psi_old = sp.copy(psi)
    omega_last = sp.copy(omega)

    # Calculo de vz e vr
    vz = sp.multiply(MinvLump, sp.dot(Gy2, psi))
    vr = -1.0 * sp.multiply(MinvLump, sp.dot(Gx2, psi))

    # vz = sp.multiply(MinvLump, sp.dot(Gy, psi))
    # vr = -1.0 * sp.multiply(MinvLump, sp.dot(Gx, psi))

# ----------------- Fim de Loop -------------------
# -------------------------------------------------
