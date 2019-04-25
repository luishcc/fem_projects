# -*- coding: utf-8 -*-
import scipy as sp
import GMesh
import libmalha
import random

# To Do:
# search1D - Identify points p out of the mesh as if they were on the closest boundary
# search2D - Identify points p out of the mesh as if they were on the closest boundary
#          - Go to next point and check neighbour elements

malha = GMesh.GMesh("test.msh")

#x, ien = libmalha.lin1d(100, 5)
x, y, ien = libmalha.lin2d(5, 2, 3, 1)
ien = libmalha.ret_tri(ien)

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

            if 1.0 >= Li >= 0.0 and 1.0 >= Lj >= 0.0:
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

    result = sp.ones(_np) * (-1)
    area_coord = sp.zeros((_np, 3))

    for i in range(_np):
        px = _xd[i]
        py = _yd[i]
        end = 0

        while end == 0:
            for e in _neighbours[i]:
                v1 = _ien[e, 0]
                v2 = _ien[e, 1]
                v3 = _ien[e, 2]

                x1 = _x[v1]
                x2 = _x[v2]
                x3 = _x[v3]

                y1 = _y[v1]
                y2 = _y[v2]
                y3 = _y[v3]

                a1 = x1 * y2 - x2 * y1
                a2 = x3 * y1 - x1 * y3
                a3 = x2 * y3 - x3 * y2
                Area_x2 = a1 + a2 + a3
                Area_x2_inv = 1.0 / Area_x2

                b1 = y2 - y3
                b2 = y3 - y1
                b3 = y1 - y2
                c1 = x3 - x2
                c2 = x1 - x3
                c3 = x2 - x1

                Li = Area_x2_inv * (a1 + b1 * px + c1 * py)
                Lj = Area_x2_inv * (a2 + b2 * px + c2 * py)
                Lk = Area_x2_inv * (a3 + b3 * px + c3 * py)

                if 1.0 >= Li >= 0.0 and 1.0 >= Lj >= 0.0 and 1.0 >= Lk >= 0.0:
                    result[i] = e
                    area_coord[i][0] = Li
                    area_coord[i][1] = Lj
                    area_coord[i][2] = Lk
                    end = 1
                    break
                else:
                    dist[0] = sp.sqrt((x1 - px) ** 2 + (y1 - py) ** 2)
                    dist[1] = sp.sqrt((x2 - px) ** 2 + (y2 - py) ** 2)
                    dist[2] = sp.sqrt((x3 - px) ** 2 + (y3 - py) ** 2)
                    end = 0

                end = 1
    return result, area_coord

dt = 0.0001
v = sp.ones((len(x),2))*0.1

for i in range(len(x)):
    v[i][0] = random.randint(-101,101) * 0.001
    v[i][1] = random.randint(-101,101) * 0.001

xd = x - v[:, 0] * dt
yd = y - v[:, 1] * dt

neib = libmalha.neighbourElements(len(x), ien)

res, area_coord = search2D(len(x), neib, ien, x, y, xd, yd)

print res, "\n", area_coord

for e in neib[1]:
    print e

exit()



