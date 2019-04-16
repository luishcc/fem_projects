import GMesh
import scipy as sp

malha = GMesh.GMesh("exemplo.msh")

#print malha.nodes_in_physical_groups[1][2]

#print malha.physical_groups_dirichlet
#print malha.nodes_in_physical_groups[1]


print malha.Boundary_Neumann
print malha.neumann_points
