# -*- coding: utf-8 -*-
import GMesh as gm
import InOut as IO
import scipy as sp
import os
import semiLagrangean as sl
cwd = os.getcwd()

#mesh_file = "test"
mesh_file = "test-s"
savename = mesh_file


mesh = gm.GMesh("mesh/" + mesh_file + ".msh")

num_nodes = len(mesh.X)
num_ele = len(mesh.IEN)

neighbour_ele, neighbour_nodes = sl.neighbourElements2(num_nodes, mesh.IEN)

print mesh.Boundary_Nodes

def smoothMesh(_neighbour_nodes, _mesh, _x, _y, _dt, delT = 1):

    xx = sp.copy(_x)
    yy = sp.copy(_y)

    for i in range(len(_neighbour_nodes)):

        flag = 0
        for k in _mesh.Boundary_Nodes:
            if i == k:
                flag = 1
                break

        if flag == 1:
            continue

        vertex_position = sp.array([xx[i], yy[i]])
        nghN = _neighbour_nodes[i]
        num_nghb = len(nghN)
        distance_vectors = sp.zeros((num_nghb, 2))
        displacement_vector = sp.zeros(2)
        for j in range(num_nghb):
            index = nghN[j]
            distance_vectors[j][0] = xx[index] - vertex_position[0]
            distance_vectors[j][1] = yy[index] - vertex_position[1]
            displacement_vector += distance_vectors[j]
        displacement_vector *= (1./num_nghb)
        displacement_velocity = displacement_vector / delT


        #new_position = vertex_position + displacement_vector
        new_position = vertex_position + displacement_velocity * _dt

        xx[i] = new_position[0]
        yy[i] = new_position[1]


      #  print nghN, "\n", distance_vectors, "\n", displacement_vector, "\n", new_position



    return xx, yy

def weightedSmoothMesh(_neighbour_nodes, _mesh, _x, _y, _dt, delT = 1):

    xx = sp.copy(_x)
    yy = sp.copy(_y)

    for i in range(len(_neighbour_nodes)):

        flag = 0
        for k in _mesh.Boundary_Nodes:
            if i == k:
                flag = 1
                break
        if flag == 1:
            continue

        vertex_position = sp.array([xx[i], yy[i]])
        nghN = _neighbour_nodes[i]
        num_nghb = len(nghN)
        neighbour_positions = sp.zeros((num_nghb, 2))
        displacement_vector = sp.zeros(2)
        sumWeight = 0
        for j in range(num_nghb):
            index = nghN[j]
            neighbour_positions[j][0] = xx[index]
            neighbour_positions[j][1] = yy[index]
            distance = sp.sqrt((xx[index] - vertex_position[0])**2 +
                               (yy[index] - vertex_position[1])**2)
#            distance = 1

            #displacement_vector += (1./distance) * neighbour_positions[j]
            #sumWeight += (1./distance)

            displacement_vector += (distance) * neighbour_positions[j]
            sumWeight += (distance)

        displacement_vector = (displacement_vector/sumWeight) - vertex_position

        displacement_velocity = displacement_vector / delT

        new_position = vertex_position + displacement_velocity * _dt

        xx[i] = new_position[0]
        yy[i] = new_position[1]


      #  print nghN, "\n", distance_vectors, "\n", displacement_vector, "\n", new_position



    return xx, yy


x2 = sp.copy(mesh.X)
y2 = sp.copy(mesh.Y)
vtk_t = IO.InOut(x2, y2, mesh.IEN, num_nodes, num_ele, None, None, None, None, None, None, None)
vtk_t.saveVTK(cwd + "/results", savename + str(0))

iter = 1000
for t in range(1, iter):

    #x2, y2 = smoothMesh(neighbour_nodes, mesh, x2, y2, 0.001, delT=0.1)
    x2, y2 = weightedSmoothMesh(neighbour_nodes, mesh, x2, y2, 0.001, delT=0.1)

    vtk_t = IO.InOut(x2, y2, mesh.IEN, num_nodes, num_ele, None, None, None, None, None, None, None)
    vtk_t.saveVTK(cwd + "/results", savename + str(t))


#print neighbour_nodes, "\n", neighbour_ele



