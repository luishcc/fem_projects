# -*- coding: utf-8 -*-
import scipy as sp
import libmalha as lm

tempo = 100
dt = 0.01
Re = 10

dtinv = 1.0/dt
reinv = 1.0/Re

    # Definição da malha

x, y, ien = libmalha.lin2d(5, 50, 1, 10)
ien = libmalha.ret_tri(ien)
np = len(x)
Nelem = len(ien)

    #Implementando Condicao inicial e velocidade

v = sp.zeros((np, 2))
U = 1.5
a = - 4 * U
b = 4 * U
for i in range(np):
    v[i][0] = a * (y[i]**2) + b * y[i]

x_d = x - v[:, 0] * dt
y_d = y - v[:, 1] * dt


T_ini = sp.zeros(np)
for i in range(np):
    T_ini[i] = 10 * (x[i]-0.5)**2 + 10 * (y[i]-0.5)**2


# Assembly

def fem_matrix(_x, _y, _numele, _numnode, _ien):
    k_local = sp.zeros((3, 3), dtype="float64")
    m_local = sp.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype="float64")
    a = sp.zeros(3, dtype="float64")
    b = sp.zeros(3, dtype="float64")
    c = sp.zeros(3, dtype="float64")
    yy = sp.zeros(3, dtype="float64")
    xx = sp.zeros(3, dtype="float64")
    k_global = sp.zeros((_numnode, _numnode), dtype="float64")
    m_global = sp.zeros((_numnode, _numnode), dtype="float64")

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

        for i_local in range(3):
            i_global = _ien[elem, i_local]
            for j_local in range(3):
                j_global = _ien[elem, j_local]
                k_global[i_global, j_global] += k_local[i_local, j_local]
                m_global[i_global, j_global] += m_local[i_local, j_local] * (Area/12.)


    return  k_global, m_global

K, M = fem_matrix(x, y, Nelem, np, ien)

kre = K * reinv
mdt = M * dtinv

T_i = sp.copy(T_ini)



exit()

# buscar elemento de x_d e y_d
# TO DO

# Time Loop

for t in range(tempo):

    # interpolar T_i para descobrir T_d
    # TO DO

    # sistema

    LHS = mdt + kre
    RHS = sp.dot((mdt) * T_d)

    #cc

    T_i = sp.linalg.solve(LHS, RHS)

