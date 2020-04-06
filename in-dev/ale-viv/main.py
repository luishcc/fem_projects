# -*- coding: utf-8 -*-
import scipy as sp
from scipy import linalg
import InOut as Io
import GMesh as Gm
import meshio
import os
import semiLagrangean as sl

cwd = os.getcwd()

arquivo = "viv"

malha = Gm.GMesh("mesh/"+arquivo+".msh")
x = malha.X
y = malha.Y
ien = malha.IEN
nodes = len(x)
num_ele = len(ien)

dt = 0.001
tempo = 200
Re = 100
v_in = Re
psi_top = v_in * 10

p_lagrange = 0.0
p_smooth = 0.
p_wave = 0.0

# ---------------------------------------
# Wz, Psi e velocidade inicial
# ---------------------------------------

Psi_new = sp.zeros(nodes, dtype="float64")
Wz_new = sp.zeros(nodes, dtype="float64")
vx = sp.zeros(nodes, dtype="float64") + v_in
vy = sp.zeros(nodes, dtype="float64")

# ---------------------------------------
# Montagem de matrizes
# ---------------------------------------

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
    gx_global = sp.zeros((_numnode, _numnode), dtype="float64")
    gy_global = sp.zeros((_numnode, _numnode), dtype="float64")

    for elem in range(_numele):
        for i in range(3):
            xx[i] = _x[_ien[elem, i]]
            yy[i] = _y[_ien[elem, i]]

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
                k_global[i_global, j_global] += k_local[i_local, j_local]
                m_global[i_global, j_global] += m_local[i_local, j_local]* (Area/12.)
                gx_global[i_global, j_global] += gx_local[i_local, j_local]
                gy_global[i_global, j_global] += gy_local[i_local, j_local]


    return  k_global, m_global, gx_global, gy_global

K, M, Gx, Gy = fem_matrix(x, y, num_ele, nodes, ien)

Boundary = sp.zeros(nodes)
lista = []

psi_bc = sp.zeros(nodes)

cylinder = Gm.get_boundary_with_tag("mesh/"+arquivo+".msh", 5)
center = 0
len_cylinder = len(cylinder)
for i in cylinder:
    center += y[i] / len_cylinder

for i in range(nodes):
    if y[i] == 0.0 :
        Boundary[i] = i
        psi_bc[i] = 0
    elif y[i] == 10.0 :
        Boundary[i] = i
        psi_bc[i] = psi_top
    elif i in cylinder:
        Boundary[i] = i
        psi_bc[i] = 0.1 * center * psi_top
    elif x[i] == 0.0 or x[i] == 32.5:
        Boundary[i] = i
    else:
        lista.append(i)

Boundary = sp.delete(Boundary, lista, axis=0)
num_bc = len(Boundary)

# ---------------------------------------
# Condições de contorno e Inicial
# ---------------------------------------

bc_omega = sp.zeros(nodes)
for i in Boundary:
    j = int(i)
    if x[j] == 0.0:
        vx[j] = v_in
    if x[j] != 32.5:
        #vx[j] = 0.0
        vy[j] = 0.
        if j in cylinder:
            vx[j] = 0




Minv = sp.linalg.inv(M)
Wz_old = sp.dot(Minv, (sp.dot(Gx, vy) - sp.dot(Gy, vx)))

K_psi = sp.copy(K)
ccpsi = sp.zeros(nodes)

for i in range(num_bc):
    index = int(Boundary[i])
    value = psi_bc[index]
    if y[index] == 0.0 or y[index] == 10.0 or (index in cylinder):
        ccpsi[index] -= value * K[j, index]
        for j in range(nodes):
            if j != index:
                K_psi[index, j] = 0
                K_psi[j, index] = 0
            else:
                K_psi[index, j] = 1

F_psi = sp.dot(M, Wz_old) + ccpsi
for i in range(num_bc):
    index = int(Boundary[i])
    F_psi[index] = psi_bc[index]

Psi_old = sp.linalg.solve(K_psi, F_psi)

sp.random.seed(1)

# ----------------------------------------------------------
# ---------------------- Loop No Tempo ------------------------

neighbour_ele, neighbour_nodes = sl.neighbourElements2(nodes, ien)

for t in range(0, tempo-1):
    print("Solving System " + str((float(t)/(tempo-1))*100) + "%")

    vx_smooth, vy_smooth = Gm.smoothMesh(neighbour_nodes, malha, x, y, dt)
    vx_wave = sp.random.rand(nodes) * 0.2*sp.cos(2*sp.pi*5*t*dt)
    vy_wave = sp.random.rand(nodes) * 0.2*sp.cos(2*sp.pi*5*t*dt)

    vxAle = p_lagrange * vx + p_smooth * vx_smooth + p_wave * vx_wave
    vyAle = p_lagrange * vy + p_smooth * vy_smooth + p_wave * vy_wave

    for i in range(num_bc):
        index = int(Boundary[i])
        vxAle[index] = 0
        vyAle[index] = 0
        #vy[index] = 0.0
        if x[index] == 0:
            vx[index] = v_in
        if x[index] != 32.5:
            #vx[index] = 0.0
            vy[index] = 0.0
            if index in cylinder:
                vx[index] = 0


    vx_sl = vx - vxAle
    vy_sl = vy - vyAle
    Wz_dep = sl.Linear2D(nodes, neighbour_ele, ien, x, y, vx_sl, vy_sl, dt, Wz_old)

    # Solução de Wz e Psi
    x = x + vxAle * dt
    y = y + vyAle * dt
    K, M, Gx, Gy = fem_matrix(x, y, num_ele, nodes, ien)

    Minv = linalg.inv(M)

    # B.C. Vorticidade
    Wcc = sp.dot(Minv, (sp.dot(Gx, vy) - sp.dot(Gy, vx)))
    ccomega = sp.zeros(nodes)


    LHS = M / dt + K / Re
    LHS_omega = sp.copy(LHS)

    for i in range(num_bc):
        index = int(Boundary[i])
        value = Wcc[index]
        for j in range(nodes):
            ccomega[j] -= value * LHS[j, index]
            if j != index:
                LHS_omega[index, j] = 0
                LHS_omega[j, index] = 0
            else:
                LHS_omega[index, j] = 1

    F_omega = sp.dot(M / dt, Wz_dep) + ccomega / Re

    for i in range(num_bc):
        index = int(Boundary[i])
        F_omega[index] = Wcc[index]

    Wz_new = sp.linalg.solve(LHS_omega, F_omega)

    K_psi = sp.copy(K)
    ccpsi = sp.zeros(nodes)

    for i in range(num_bc):
        index = int(Boundary[i])
        value = psi_bc[index]
        if y[index] == 0.0 or y[index] == 10.0 or (index in cylinder):
            for j in range(nodes):
                ccpsi[j] -= value * K[j, index]
                if j != index:
                    K_psi[index, j] = 0
                    K_psi[j, index] = 0
                else:
                    K_psi[index, j] = 1

    F_psi = sp.dot(M, Wz_new) + ccpsi
    for i in range(num_bc):
        index = int(Boundary[i])
        F_psi[index] = psi_bc[index]

    Psi_new = sp.linalg.solve(K_psi, F_psi)

    # Salvar VTK
    vtk = Io.InOut(x, y, ien, len(x), len(ien), Psi_old, Wz_old, Wz_dep,
                    None, None, vx, vy, vx_sl, vy_sl)
    vtk.saveVTK(cwd+"/results", arquivo + str(t+1))

    Psi_old = sp.copy(Psi_new)
    Wz_old = sp.copy(Wz_new)

    # Calculo de Vx e Vy
    vx = sp.dot(Minv, sp.dot(Gy, Psi_new))
    vy = -1.0 * sp.dot(Minv, sp.dot(Gx, Psi_new))

#----------------- Fim de Loop -------------------
#-------------------------------------------------
