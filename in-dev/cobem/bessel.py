# -*- coding: utf-8 -*-
import GMesh as gm
import InOut as IO
import scipy as sp
from scipy import special
import os


cwd = os.getcwd()

# -------------------------------------------------------
#     Reading Mesh
# -------------------------------------------------------

mesh_file = "poiseuille"
savename = "bessel"
fluid_mesh = gm.GMesh("mesh/" + mesh_file + "-fld.msh")

x_fluid = fluid_mesh.X
y_fluid = fluid_mesh.Y
ien_fluid = fluid_mesh.IEN
nodes_fluid = len(x_fluid)
num_ele_fluid = len(ien_fluid)

temp = sp.zeros(nodes_fluid)

maxsum = 100000
for i in range(nodes_fluid):
	print i, " / ", nodes_fluid
	z = x_fluid[i]
	r = y_fluid[i]
	i0 = sp.special.iv(2, 1)
	ir =  sp.special.iv(2, r)
	for n in range(1,maxsum):
		ln = n*(sp.pi/5.)	
		c1 = (2./ln) * (sp.sin(2.5*ln))**2
		c2 = 2.5 - sp.sin(10*ln)/(4*ln)
		cn = c1 / (i0 * c2)
		temp[i] += cn * ir * sp.sin(z*ln) 
	

# Salvar VTK
vtk_t = IO.InOut(x_fluid, y_fluid, ien_fluid, nodes_fluid, num_ele_fluid, temp, None, None
                 , None, None, None, None)
vtk_t.saveVTK(cwd+"/resultsTest", savename + "Test")

