# -*- coding: utf-8 -*-
import scipy as sp
import GMesh
import libmalha
import random

# To Do:
# search1D - Identify points p out of the mesh as if they were on the closest boundary


malha = GMesh.GMesh("test.msh")

x, ien = libmalha.lin1d(100, 5)

def search1D(_np, _ne, _ien, _x, _xd, order='1'):
    print "searching element 1D with order = ", order

    result = sp.zeros(_np)
    area_coord = sp.zeros((_np,2))


    for i in range(0,_np):
        p = _xd[i]

        for e in range(0, _ne):
            x1 = _x[_ien[e, 0]-1]
            x2 = _x[_ien[e, 1]-1]

            l1 = abs(x2 - p)
            l2 = abs(x1 - p)
            l = abs(x2 - x1)
            linv = 1.0 / l

            Li = l1 * linv
            Lj = l2 * linv

            if 1.0 >= Li >= 0 and 1.0 >= Lj >= 0:
                result[i] = e
                area_coord[i][0] = Li
                area_coord[i][1] = Lj

    return result, area_coord
#
# dt = 0.021
# v = sp.ones(len(x))
# for i in range(len(x)):
#     v[i] = random.randint(-101,101) * 0.01
#
# xd = x - v * dt
#
# res, area_coord = search1D(len(x), len(ien), ien, x, xd)
#
# print res, area_coord


def search2D(_np, _neighbours, _ien, _x, _y, _xd, _yd, order='1'):
    result = sp.zeros(_np)
    area_coord = sp.zeros(_np)

    for i in range(_np):
        px = _xd[i]
        py = _yd[i]

        for e in _neighbours:
            v1 = _ien[e, 0]
            v2 = _ien[e, 1]
            v3 = _ien[e, 2]

            x1 = _x[v1]
            x2 = _x[v2]
            x3 = _x[v3]

            y1 = _y[v1]
            y2 = _y[v2]
            y3 = _y[v3]


    return