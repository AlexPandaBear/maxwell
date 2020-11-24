print("----------------------------------------------------------------------")
print("---------------------------- MAXWELL TEST ----------------------------")
print("----------------------------------------------------------------------\n")




print("--- [ LOADING MODULES ] ----------------------------------------------")

print("Matplotlib", end=" ")

import matplotlib.pyplot as plt
import matplotlib.animation as anim
import matplotlib.cm as cm
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D

print("[ OK ]")

print("NumPy", end=" ")

import numpy as np
from itertools import combinations

print("[ OK ]")

print("Maxwell", end=" ")

mxll_path = "../build/install/lib/python"
import sys

if not (mxll_path in sys.path):
	sys.path.append(mxll_path)

import _maxwell as mxll

print("[ OK ]")

print("--------------------------------------------------------- [ DONE ] ---\n")




print("--- [ LOADING PARAMETERS ] -------------------------------------------")

S = mxll.ESM()


epsilon0 = 8.854187e-12
accuracy = 0.00001
max_nb_iterations = 100

S.setEpsilon0(epsilon0)
S.setAccuracy(accuracy)
S.setMaxNbIterations(max_nb_iterations)


x_min = -1.
x_max = 1.
y_min = -1.
y_max = 1.
z_min = -1.
z_max = 1.

nx = 10
ny = 10
nz = 10

S.generateCubeMesh(x_min, x_max, nx, y_min, y_max, ny, z_min, z_max, nz)


tot_nb_nodes = nx*ny*nz

rho = mxll.ScalarField(tot_nb_nodes, 0.)
eps_r = mxll.ScalarField(tot_nb_nodes, 1.)

for i in range(tot_nb_nodes):
	xyz = S.getNodeXYZ(i)
	x = xyz.getX()
	y = xyz.getY()
	z = xyz.getZ()

	if x**2 + y**2 + z**2 > 4.:
		rho.setValue(i, 1.)

S.defineRhoField(rho)
S.defineEpsilonRField(eps_r)

print("--------------------------------------------------------- [ DONE ] ---\n")




print("--- [ COMPUTING ] ----------------------------------------------------")

S.simulate()

print("--------------------------------------------------------- [ DONE ] ---\n")




print("--- [ EXTRACTING DATA ] ----------------------------------------------")

print("Extracting discretisation...")

nb_nodes = S.getNbNodes()
M = np.asarray(S.getMesh())

print("Extracting fields...")

E = np.asarray(S.getE())

print("--------------------------------------------------------- [ DONE ] ---\n")




print("--- [ PLOTTING ] -----------------------------------------------------")

print("Plotting electromagnetic field...")
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.quiver(M[:,0], M[:,1], M[:,2], E[:,0], E[:,1], E[:,2], normalize=True)
#fig.colorbar(cm.ScalarMappable(norm=None, cmap='jet'), ax=ax)
plt.show()

print("--------------------------------------------------------- [ DONE ] ---\n")



print("----------------------------------------------------------------------")
print("----------------------------------------------------------------------")
print("----------------------------------------------------------------------")