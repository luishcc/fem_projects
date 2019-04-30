# -*- coding: utf-8 -*-
import scipy as sp
import GMesh
import libmalha
from scipy import linalg
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
    outside_point = []

    for i in range(_np):
        px = _xd[i]
        py = _yd[i]
        end = 0
        j = i
        dist = []

        while end == 0:
            for e in _neighbours[j]:

                print e , i
                x1 = _x[_ien[e, 0]]
                x2 = _x[_ien[e, 1]]
                x3 = _x[_ien[e, 2]]

                y1 = _y[_ien[e, 0]]
                y2 = _y[_ien[e, 1]]
                y3 = _y[_ien[e, 2]]

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
                    dist.append([_ien[e, 0], sp.sqrt((x1 - px) ** 2 + (y1 - py) ** 2)])
                    dist.append([_ien[e, 1], sp.sqrt((x2 - px) ** 2 + (y2 - py) ** 2)])
                    dist.append([_ien[e, 2], sp.sqrt((x3 - px) ** 2 + (y3 - py) ** 2)])
                    end = 0

            if end == 1:
                break
            else:
                min_dist = min(dist, key=lambda k: k[1])
                last_node = j
                j = min_dist[0]

                if last_node == j and end == 0:
                    outside_point.append([i, last_node])
                    end = 1
                    break

    return result, area_coord


#-------------------------------
#-------------------------------
#-------------------------------
def Linear2D(_npoints, _neighbors_elements, _IEN, _xn, _yn, _xd, _yd):
    xd = _xd
    yd = _yd

    result = sp.ones(_npoints) * (-1)
    area_coord = sp.zeros((_npoints, 3))


    for i in range(0, _npoints):
        x = float(xd[i])
        y = float(yd[i])

        node = i
        length = []
        breaking = 0

        while breaking == 0:
            for e in _neighbors_elements[node]:
                v1 = _IEN[e][0]
                v2 = _IEN[e][1]
                v3 = _IEN[e][2]

                x1 = float(_xn[v1])
                x2 = float(_xn[v2])
                x3 = float(_xn[v3])

                y1 = float(_yn[v1])
                y2 = float(_yn[v2])
                y3 = float(_yn[v3])

                A = sp.array([[x1, x2, x3],
                              [y1, y2, y3],
                              [1.0, 1.0, 1.0]])

                b = sp.array([x, y, 1.0])

                alpha = sp.linalg.solve(A, b)

                if sp.all(alpha >= 0.0) and sp.all(alpha <= 1.0):

                    A1 = 0.5 * sp.linalg.det(sp.array([[1, x, y],
                                                       [1, x2, y2],
                                                       [1, x3, y3]]))

                    A2 = 0.5 * sp.linalg.det(sp.array([[1, x1, y1],
                                                       [1, x, y],
                                                       [1, x3, y3]]))

                    A3 = 0.5 * sp.linalg.det(sp.array([[1, x1, y1],
                                                       [1, x2, y2],
                                                       [1, x, y]]))

                    At = 0.5 * sp.linalg.det(sp.array([[1, x1, y1],
                                                       [1, x2, y2],
                                                       [1, x3, y3]]))

                    Li = A1 / At
                    Lj = A2 / At
                    Lk = A3 / At

                    area_coord[i][0] = Li
                    area_coord[i][1] = Lj
                    area_coord[i][2] = Lk
                    result[i] = e

                    breaking = 1
                    break


                else:
                    x_a = x1 - x
                    x_b = x2 - x
                    x_c = x3 - x

                    y_a = y1 - y
                    y_b = y2 - y
                    y_c = y3 - y

                    length1 = sp.sqrt(x_a ** 2 + y_a ** 2)
                    length2 = sp.sqrt(x_b ** 2 + y_b ** 2)
                    length3 = sp.sqrt(x_c ** 2 + y_c ** 2)

                    a_1 = [v1, length1]
                    a_2 = [v2, length2]
                    a_3 = [v3, length3]

                    length.append(a_1)
                    length.append(a_2)
                    length.append(a_3)

                    breaking = 0

            if breaking == 1:
                break

            else:
                length_min = min(length, key=lambda k: k[1])
                node1 = node
                node = length_min[0]

                # outside domain
                if node == node1 and breaking == 0:

                    breaking = 1
                    break

    return result, area_coord
#-------------------------------




dt = 0.0001
v = sp.ones((len(x),2))*0.1
#
# for i in range(len(x)):
#     v[i][0] = random.randint(-101,101) * 0.001
#     v[i][1] = random.randint(-101,101) * 0.001

xd = x - v[:, 0] * dt
yd = y - v[:, 1] * dt

neib = libmalha.neighbourElements(len(x), ien)

res, area_coord = search2D(len(x), neib, ien, x, y, xd, yd)
res2, area2 = Linear2D(len(x), neib, ien, x, y, xd, yd)

print res, "\n", area_coord

print res2, "\n", area2

#
# for e in neib[1]:
#     print e

exit()



