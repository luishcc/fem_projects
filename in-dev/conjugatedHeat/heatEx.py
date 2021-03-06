# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import InOut as Io
import GMesh as Gm
import os

cwd = os.getcwd()


# -------------------------------------------------------
#     Reading Mesh
# -------------------------------------------------------

arquivo = "heat"

malha_total = Gm.GMesh(arquivo+".msh")

x_total = malha_total.X
y_total = malha_total.Y
ien_total = malha_total.IEN
nodes_total = len(x_total)
num_ele_total = len(ien_total)

malha_fldH = Gm.GMesh(arquivo+"-fldh.msh")

x_fldH = malha_fldH.X
y_fldH = malha_fldH.Y
ien_fldH = malha_fldH.IEN
nodes_fldH = len(x_fldH)
num_ele_fldH = len(ien_fldH)

malha_fldC = Gm.GMesh(arquivo+"-fldc.msh")

x_fldC = malha_fldC.X
y_fldC = malha_fldC.Y
ien_fldC = malha_fldC.IEN
nodes_fldC = len(x_fldC)
num_ele_fldC = len(ien_fldC)



for i in range(nodes_total):
    if y_total[i] > 0.3:
        print y_total[i]
    if y_total[i] < 0.0:
        print y_total[i]

for i in range(nodes_fldH):
    if y_fldH[i] > 0.3:
        print y_fldH[i], "  -H"
    if y_fldH[i] < 0.2:
        print y_fldH[i], "  -H"

for i in range(len(malha_fldH.dirichlet_points)):
    index = int(malha_fldH.dirichlet_points[i][0]-1)
    if y_fldH[index]!= 0.2 and y_fldH[index] != 0.3:
        print y_fldH[index], " - bound H"

for i in range(nodes_fldC):
    if y_fldC[i] > 0.1:
        print y_fldC[i], "  -C"
    if y_fldC[i] < 0.0:
        print y_fldC[i], "  -C"

for i in range(len(malha_fldC.dirichlet_points)):
    index = int(malha_fldC.dirichlet_points[i][0]-1)
    if y_fldC[index]!= 0.0 and y_fldC[index] != 0.1:
        print y_fldC[index], " - bound C"


exit()



Vel_points_convertH = sp.zeros(nodes_total) -1
for i in range(nodes_total):
    for j in range(nodes_fldH):
        if (x_total[i] == x_fldH[j]) and (y_total[i] == y_fldH[j]):
            Vel_points_convertH[i] = j

Vel_points_convertC = sp.zeros(nodes_total) - 1
for i in range(nodes_total):
    for j in range(nodes_fldC):
        if (x_total[i] == x_fldC[j]) and (y_total[i] == y_fldC[j]):
            Vel_points_convertC[i] = j



# -------------------------------------------------------
#     Defining Parameters
# -------------------------------------------------------

dt = 0.01
tempo = 700

Ratio_k = 5

Ni_H = 0.00015
alfa_fldH = 0.0002

Ni_C = Ni_H
alfa_fldC = alfa_fldH

alfa_solid = Ratio_k * alfa_fldH


# ---------------------------------------
# Wz, Psi, T e velocidade inicial
# ---------------------------------------

Psi_new_h = sp.zeros(nodes_fldH, dtype="float64")
Wz_new_h = sp.zeros(nodes_fldH, dtype="float64")
vx_h = sp.zeros(nodes_fldH, dtype="float64")
vy_h = sp.zeros(nodes_fldH, dtype="float64")

Psi_new_c = sp.zeros(nodes_fldC, dtype="float64")
Wz_new_c = sp.zeros(nodes_fldC, dtype="float64")
vx_c = sp.zeros(nodes_fldC, dtype="float64")
vy_c = sp.zeros(nodes_fldC, dtype="float64")


vx_total = sp.zeros(nodes_total, dtype="float64")
vy_total = sp.zeros(nodes_total, dtype="float64")

Temp_old = sp.ones(nodes_total, dtype="float64")*550
for i in range(len(malha_total.dirichlet_points)):
    index = int(malha_total.dirichlet_points[i][0]-1)
    value = malha_total.dirichlet_points[i][1]
    Temp_old[index] = value


# ---------------------------------------
# Montagem de matrizes
# ---------------------------------------

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

# --------------------------------------
# Total Matrices (Heat equation)


Alfa_total = sp.zeros(nodes_total)
for i in range(nodes_total):
    if y_total[i] < 0.1:
        Alfa_total[i] = alfa_fldH
    if y_total[i] > 0.2:
        Alfa_total[i] = alfa_fldH
    else:
        Alfa_total[i] = alfa_solid

K, M, Gx, Gy = fem_matrix(x_total, y_total, num_ele_total, nodes_total, ien_total, Alfa_total)

Mdt = M/dt
# --------------------------------------


# --------------------------------------
# Fluid Region (Psi, Omega)

# Hot Fluid
Alfa_fluid_h = sp.ones(nodes_fldH)

K_fld_h, M_fld_h, Gx_fld_h, Gy_fld_h = fem_matrix(x_fldH, y_fldH, num_ele_fldH, nodes_fldH, ien_fldH, Alfa_fluid_h)

K_fld_ni_h = K_fld_h * Ni_H
Mdt_fld_h = M_fld_h/dt

MLump_h = sp.zeros((nodes_fldH, nodes_fldH))
for i in range(nodes_fldH):
    for j in range(nodes_fldH):
        MLump_h[i, i] += M_fld_h[i, j]
MinvLump_h = linalg.inv(MLump_h)

K_psi_h = sp.copy(K_fld_h)
ccpsi_h = sp.zeros(nodes_fldH)

# Cold Fluid
Alfa_fluid_c = sp.ones(nodes_fldC)

K_fld_c, M_fld_c, Gx_fld_c, Gy_fld_c = fem_matrix(x_fldC, y_fldC, num_ele_fldC, nodes_fldC, ien_fldC, Alfa_fluid_c)

K_fld_ni_c = K_fld_c * Ni_C
Mdt_fld_c = M_fld_c/dt

MLump_c = sp.zeros((nodes_fldC, nodes_fldC))
for i in range(nodes_fldC):
    for j in range(nodes_fldC):
        MLump_c[i, i] += M_fld_c[i, j]
MinvLump_c = linalg.inv(MLump_c)

K_psi_c = sp.copy(K_fld_c)
ccpsi_c = sp.zeros(nodes_fldC)
# --------------------------------------

# ---------------------------------------
# Condições de contorno e Inicial
# ---------------------------------------

dirichlet_len_total = len(malha_total.dirichlet_points)
dirichlet_len_fluid_h = len(malha_fldH.dirichlet_points)
dirichlet_len_fluid_c = len(malha_fldC.dirichlet_points)

# Hot
Fluid_Boundary_h = sp.zeros(nodes_fldH)
lista = []
for i in range(nodes_fldH):
    if x_fldH[i] == 0.0 or y_fldH[i] == 0.2 or y_fldH[i] == 0.3 or x_fldH[i] == 1.0:
        Fluid_Boundary_h[i] = i
    else:
        lista.append(i)
Fluid_Boundary_h = sp.delete(Fluid_Boundary_h, lista, axis=0)

num_omega_h = len(Fluid_Boundary_h)

#Cold
Fluid_Boundary_c = sp.zeros(nodes_fldC)
lista1 = []
for i in range(nodes_fldC):
    if x_fldC[i] == 0.0 or y_fldC[i] == 0.0 or y_fldC[i] == 0.1 or x_fldC[i] == 1.0:
        Fluid_Boundary_c[i] = i
    else:
        lista1.append(i)
Fluid_Boundary_c = sp.delete(Fluid_Boundary_c, lista1, axis=0)

num_omega_c = len(Fluid_Boundary_c)


# --------------------------------------
# Psi K matrix with Dirichlet BC -- K_psi

# Hot
for i in range(dirichlet_len_fluid_h):
    index = int(malha_fldH.dirichlet_points[i][0] - 1)
    value = malha_fldH.dirichlet_points[i][1]
    for j in range(nodes_fldH):
        ccpsi_h[j] -= value * K_fld_h[j, index]
        if j != index:
            K_psi_h[index, j] = 0
            K_psi_h[j, index] = 0
        else:
            K_psi_h[index, j] = 1

# Cold
for i in range(dirichlet_len_fluid_c):
    index = int(malha_fldC.dirichlet_points[i][0] - 1)
    value = malha_fldC.dirichlet_points[i][1]
    for j in range(nodes_fldC):
        ccpsi_c[j] -= value * K_fld_c[j, index]
        if j != index:
            K_psi_c[index, j] = 0
            K_psi_c[j, index] = 0
        else:
            K_psi_c[index, j] = 1

# --------------------------------------


# --------------------------------------
# Velocity

# Hot
for i in Fluid_Boundary_h:
    j = int(i)
    vy_h[j] = 0.0
    if x_fldH[j] == 0.0:
        vx_h[j] = 0.2
    if y_fldH[j] == 0.2 or y_fldH[j] == 0.3:
        vx_h[j] = 0.0

# Cold
for i in Fluid_Boundary_c:
    j = int(i)
    vy_c[j] = 0.0
    if x_fldC[j] == 1.0:
        vx_c[j] = -0.1
    if y_fldC[j] == 0.1 or y_fldC[j] == 0.0:
        vx_c[j] = 0.0
# --------------------------------------


# Hot
Wz_old_h = sp.dot(MinvLump_h, (sp.dot(Gx_fld_h, vy_h) - sp.dot(Gy_fld_h, vx_h)))

F_psi_h = sp.dot(M_fld_h, Wz_old_h) + ccpsi_h
for i in range(dirichlet_len_fluid_h):
    index = int(malha_fldH.dirichlet_points[i][0] - 1)
    F_psi_h[index] = malha_fldH.dirichlet_points[i][1]

Psi_old_h = sp.linalg.solve(K_psi_h, F_psi_h)


# Cold
Wz_old_c = sp.dot(MinvLump_c, (sp.dot(Gx_fld_c, vy_c) - sp.dot(Gy_fld_c, vx_c)))

F_psi_c = sp.dot(M_fld_c, Wz_old_c) + ccpsi_c
for i in range(dirichlet_len_fluid_c):
    index = int(malha_fldC.dirichlet_points[i][0] - 1)
    F_psi_c[index] = malha_fldC.dirichlet_points[i][1]

Psi_old_c = sp.linalg.solve(K_psi_c, F_psi_c)


T_a = sp.zeros(nodes_total)
u_a = sp.zeros(nodes_total)
Tan = 0.95
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

    # Hot
    for i in Fluid_Boundary_h:
        j = int(i)
        vy_h[j] = 0.0
        if x_fldH[j] == 0.0:
            vx_h[j] = 0.2
        if y_fldH[j] == 0.2 or y_fldH[j] == 0.3:
            vx_h[j] = 0.0

    # Cold
    for i in Fluid_Boundary_c:
        j = int(i)
        vy_c[j] = 0.0
        if x_fldC[j] == 1.0:
            vx_c[j] = -0.1
        if y_fldC[j] == 0.1 or y_fldC[j] == 0.0:
            vx_c[j] = 0.0


    for i in range(nodes_total):
        if Vel_points_convertH[i] >= 0:
            vx_total[i] = vx_h[int(Vel_points_convertH[i])]
            vy_total[i] = vy_h[int(Vel_points_convertH[i])]

    for i in range(nodes_total):
        if Vel_points_convertC[i] >= 0:
            vx_total[i] = vx_c[int(Vel_points_convertC[i])]
            vy_total[i] = vy_c[int(Vel_points_convertC[i])]
        # else:
        #     vx_total[i] = 0
        #     vy_total[i] = 0

    # B.C. Vorticidade Hot
    Wcc_h = sp.dot(MinvLump_h, (sp.dot(Gx_fld_h, vy_h) - sp.dot(Gy_fld_h, vx_h)))
    ccomega_h = sp.zeros(nodes_fldH)

    # Solução de Wz e Psi
    Conv_h = vx_h * Gx_fld_h + vy_h * Gy_fld_h
    LHS_Ni_h = Mdt_fld_h + K_fld_ni_h + Conv_h
    LHS_omega_h = sp.copy(LHS_Ni_h)

    for i in range(num_omega_h):
        index = int(Fluid_Boundary_h[i])
        value = Wcc_h[index]
        for j in range(nodes_fldH):
            ccomega_h[j] -= value * LHS_Ni_h[j, index]
            if j != index:
                LHS_omega_h[index, j] = 0
                LHS_omega_h[j, index] = 0
            else:
                LHS_omega_h[index, j] = 1

    F_omega_h = sp.dot(Mdt_fld_h, Wz_old_h) + ccomega_h  # - sp.dot(Conv_h, Wz_old_h)

    for i in range(num_omega_h):
        index = int(Fluid_Boundary_h[i])
        F_omega_h[index] = Wcc_h[index]

    Wz_new_h = sp.linalg.solve(LHS_omega_h, F_omega_h)

    F_psi_h = sp.dot(M_fld_h, Wz_new_h) + ccpsi_h
    for i in range(dirichlet_len_fluid_h):
        index = int(malha_fldH.dirichlet_points[i][0]-1)
        F_psi_h[index] = malha_fldH.dirichlet_points[i][1]

    Psi_new_h = sp.linalg.solve(K_psi_h, F_psi_h)

    # B.C. Vorticidade Cold
    Wcc_c = sp.dot(MinvLump_c, (sp.dot(Gx_fld_c, vy_c) - sp.dot(Gy_fld_c, vx_c)))
    ccomega_c = sp.zeros(nodes_fldC)

    # Solução de Wz e Psi
    Conv_c = vx_c * Gx_fld_c + vy_c * Gy_fld_c
    LHS_Ni_c = Mdt_fld_c + K_fld_ni_c + Conv_c
    LHS_omega_c = sp.copy(LHS_Ni_c)

    for i in range(num_omega_c):
        index = int(Fluid_Boundary_c[i])
        value = Wcc_c[index]
        for j in range(nodes_fldC):
            ccomega_c[j] -= value * LHS_Ni_c[j, index]
            if j != index:
                LHS_omega_c[index, j] = 0
                LHS_omega_c[j, index] = 0
            else:
                LHS_omega_c[index, j] = 1

    F_omega_c = sp.dot(Mdt_fld_c, Wz_old_c) + ccomega_c  # - sp.dot(Conv_c, Wz_old_c)

    for i in range(num_omega_c):
        index = int(Fluid_Boundary_c[i])
        F_omega_c[index] = Wcc_c[index]

    Wz_new_c = sp.linalg.solve(LHS_omega_c, F_omega_c)

    F_psi_c = sp.dot(M_fld_c, Wz_new_c) + ccpsi_c
    for i in range(dirichlet_len_fluid_c):
        index = int(malha_fldC.dirichlet_points[i][0]-1)
        F_psi_c[index] = malha_fldC.dirichlet_points[i][1]

    Psi_new_c = sp.linalg.solve(K_psi_c, F_psi_c)


    # Temperature
    # Conv_total = sp.dot(sp.diag(vx_total), Gx) + sp.dot(sp.diag(vy_total), Gy)
    Conv_total = vx_total * Gx + vy_total * Gy


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
    vtk_t = Io.InOut(x_total, y_total, ien_total, nodes_total, num_ele_total, Temp_old, None, None, None, None, vx_total, vy_total)
    vtk_h = Io.InOut(x_fldH, y_fldH, ien_fldH, nodes_fldH, num_ele_fldH, Psi_old_h, Wz_old_h, None, None, None, vx_h, vy_h)
    vtk_c = Io.InOut(x_fldC, y_fldC, ien_fldC, nodes_fldC, num_ele_fldC, Psi_old_c, Wz_old_c, None, None, None, vx_c, vy_c)
    vtk_t.saveVTK(cwd+"/results", arquivo + str(t+1))
    vtk_h.saveVTK(cwd+"/results", "hot" + str(t+1))
    vtk_c.saveVTK(cwd+"/results", "cold" + str(t+1))

    Psi_old_h = sp.copy(Psi_new_h)
    Wz_old_h = sp.copy(Wz_new_h)
    Psi_old_c = sp.copy(Psi_new_c)
    Wz_old_c = sp.copy(Wz_new_c)
    Temp_old = sp.copy(Temp_new)

    # Calculo de Vx e Vy
    vx_h = sp.dot(MinvLump_h, sp.dot(Gy_fld_h, Psi_new_h))
    vy_h = -1.0 * sp.dot(MinvLump_h, sp.dot(Gx_fld_h, Psi_new_h))

    vx_c = sp.dot(MinvLump_c, sp.dot(Gy_fld_c, Psi_new_c))
    vy_c = -1.0 * sp.dot(MinvLump_c, sp.dot(Gx_fld_c, Psi_new_c))

# ----------------- Fim de Loop -------------------
# -------------------------------------------------
